"""
@author: Brett Settle
@Department: UCI Neurobiology and Behavioral Science
@Lab: Parker Lab
@Date: August 6, 2015
"""
import os,sys
from BioDocks import *
from Analysis import *
from pyqtgraph.Qt import QtCore, QtGui
app = QtGui.QApplication([])

currentThread = None
currentGroup = 0

K = {'flashExists': False}

settings = \
"crop_frames=[]                         When testing parameters, crop the movie to speed up analysis.\n" \
"flashExists=true                       True if the movie contains an IR flash.\n"\
"puffDetectionMethod='slope'            ['slope' | 'simpleThreshold']\n"\
"sigmaXY=3                              filter strength.  This is sigma (in pixels) in the gaussian which will smooth the data. (default:3)\n"\
"tfiltN=[1,.015,.8]                      the order and the normalized cuttoff frequency of the temporal filter (default:[1,.015,1])\n"\
"threshold_constant=10                  decreasing this number means more puffs (and noise) will be detected (default: 8)\n"\
"min_rise_time=5                        the minimum spike rise time in frames. Must be 1 or greater (default: 2)\n"\
"max_rise_time=20                       the maximum spike rise time in frames. (default: 10)\n"\
"dilateXY=10                             how far apart puffing pixels need to be to be considered distinct puffs (in pixels) (default:2)\n"\
"dilateT=10                             (in frames) (default:10)\n"\
"minPuffSize=100                         The smallest number of puffing voxels to be considered a puff. (default:10)\n"\
"                                       for only high pass set last number to 1.  for bandpass, set it to a number less than 1.\n"\
"constants for 'analyzePuff'\n"\
"padding=30                             how many pixels do you want to pad the ROI by when fitting with gaussian? (default:40)\n"\
"wfiltN=[1,.5]                          the order and the normalized cuttoff frequency of the weak low pass filter (default: [1,.5])\n"\
"sfiltN=[1,.1]                          the order and the normalized cuttoff frequency of the strong low pass filter (default: [1,.2])\n"\
"max_fall_time=100                       this is as far after the peak as we look at. (default: 20)\n"\
"RotatedFit=true                        Set this to true to fit to 2d rotating gaussian. (default: true)\n"\
"maxSigmaForGaussianFit=20              During gaussian fit, sigma can't exceed this number. (default: 20)"

win = DockWindow(addMenu=False)
win.resize(1400, 800)
traceWidget = PlotWidget(name='Average Trace')
averageLine = pg.PlotDataItem()
traceWidget.addItem(averageLine, name='Average Line')
puffTraceWidget = PlotWidget(name='Puff Trace')
imageWidget = VideoWidget(name='Video Dock', view=ROIViewBox(roiMenu=True))
opsWidget = QWidget()
fileOps = OptionsWidget('File Properties for', [{'key': 'blacklevel', 'name': 'Blacklevel', 'value': 0}, {'key': 'crop_frames', 'name': 'Crop Frame', 'value': (0, 0)}, {'key': 'sigmaXY', 'name': 'Sigma XY', 'value': 3}, \
	{'key': 'min_rise_time', 'name': 'Minimum Rise Time', 'value': 5, 'int': True, 'step': 1}, {'key': 'dilateT', 'name': 'Frame Dilation', 'value': 10, 'int': True, 'step': 1}, {'key': 'puffDetectionMethod', 'name': 'Detection Method', 'value': ['slope', 'threshold']}, \
	{'key': 'tfiltN', 'name': 'T Filt N', 'value': [1, .015, 1]}, {'key': 'max_rise_time', 'name': 'Maximum Rise Time', 'value': 20}, {'key': 'minPuffSize', 'name': 'Minimum Puff Size', 'value': 15}, \
	{'key': 'flashExists', 'name': 'Flash Exists', 'value': True}, {'key': 'threshold_constant', 'name': 'Threshold Constant', 'value': 8}, {'key': 'dilateXY', 'name': 'XY Dilation', 'value': 5, 'int': True, 'step': 1}], shape=(3, 8))
analysisOps = OptionsWidget('Analysis Properties', [{'key': 'padding', 'name': 'Padding', 'value': 40, 'int': True, 'step': 1}, {'key': 'sfiltN', 'name': 'SFilt N', 'value': [1, .2]}, \
	{'key': 'RotatedFit', 'name': 'Rotated Fit', 'value': True}, {'key': 'wfiltN', 'name': 'W Filt N', 'value': [1, .5]}, {'key': 'max_fall_time', 'name': 'Maximum Fall Time', 'value': 100, 'int': True, 'step': 1}, \
	{'key': 'maxSigmaForGaussianFit', 'name': 'Maximum Sigma For Gaussian Fit', 'value': 20}], shape=(2, 6))


def make_opts():
	doneButton = QPushButton('Subtract Blacklevel')
	fileOps.update(K)
	analysisOps.update(K)
	doneButton.pressed.connect(ready)
	layout = QGridLayout(opsWidget)
	layout.addWidget(fileOps, 0, 0, 3, 4)
	layout.addWidget(analysisOps, 3, 0, 2, 4)
	layout.addWidget(doneButton, 6, 3)
	opsWidget.setLayout(layout)

def ready():
	global K, fileOps, analysisOps
	K.update(fileOps.getOptions())
	K.update(analysisOps.getOptions())
	K['image'] = K['image'][K['crop_frames'][0]:K['crop_frames'][1]]
	remove_flash()
	subtract_blacklevel()
	win.moveDock(puffDock, 'above', opsDock)
	beginAction.setEnabled(True)

def remove_flash():
	global K
	K['flash_frames'] = []
	if K['flashExists']:
		try:
			s, e = locate_flash(K['image'])
			K['flash_frames'] = s, e
			K['image'][s:e] = np.average(K['image'][s-10:s])
		except Exception as e:
			print('No flash located. %s' % e)
	imageWidget.setImage(K['image'])

def subtract_blacklevel():
	global K
	for i in np.arange(len(K['image'])):
		K['image'][i]=np.subtract(K['image'][i], K['blacklevel'])
		K['image'][i][K['image'][i] < 0] = 0
	imageWidget.setImage(K['image'])
	win.statusBar().showMessage("Blackleve subtracted. To continue the analysis process, press Analyze > Begin in the menubar.")

def imageChanged():
	im = imageWidget.getProcessedImage()
	global K, fileOps
	averageLine.setData(np.average(im, (2, 1)))
	K['image'] = im
	fileOps.update({'blacklevel': np.min(np.average(K['image'], 0)), 'crop_frames': [0, len(K['image'])]})

def plot_traces():
	global currentGroup
	y, x = imageWidget.location.getPosition()
	h, w = np.shape(K['image'][0])
	if not QRectF(0, 0, h, w).contains(y, x):
		return
	for i in range(len(traceWidget.getViewBox().addedItems))[::-1]:
		it = traceWidget.getViewBox().addedItems[i]
		if it.__name__.startswith(('B', 'E')):
			traceWidget.removeItem(it)
	x, y = int(x), int(y)
	puffGroup = None
	puffOriginBool=np.zeros(np.size(K['ratioed'],0)) #Falses, array of frame length
	if 'puffGroups' in K:
		puff_id, group_id=getPuffIndices(y, x)
		if group_id > 0:
			puffGroup=K['puffGroups'][group_id]
			traceWidget.setTitle("Group %i: (%.4f, %.4f)" % (group_id + 1, puffGroup.x, puffGroup.y))
			currentGroup=group_id
			for i in range(len(puffGroup.puff_ids)):
			   puff_ids=puffGroup.puff_ids[i]
			   puff=K['puffs'][puff_ids]
			   puffOriginBool[puff.t_peak]=1
		else:
			traceWidget.setTitle("(%.4f, %.4f)" % (y,x))
	else:
		traceWidget.setTitle("(%.4f, %.4f)" % (y, x))


	r=range(np.size(K['ratioed'],0))
	b=np.squeeze(K['ratioed'][:, y, x])
	e=np.squeeze(K['ratioed2'][:, y, x])
	e -= np.max(e)
	c=.4+.1*puffOriginBool
	d=.55+.1*np.squeeze(K['puffBoolean'][:, y, x])
	plotLine('b', y=b[r], pen=QPen(Qt.white))
	plotLine('e', y=e[r], pen=QPen(Qt.white))
	plotLine('c', y=c[r], pen=pg.mkPen(Qt.white, 3))
	plotLine('d', y=d[r], pen=QPen(Qt.yellow))

	if puffGroup != None: #%plot where puffing
		for i in range(len(puffGroup.puff_ids)):
			puff_ids=puffGroup.puff_ids[i]
			p=K['puffs'][puff_ids]
			t=range(int(p.t_start), int(p.t_end+1))
			plotLine("B%i" % i, x=t, y=b[t], pen=QPen(Qt.red))
			plotLine("E%i" % i, x=t, y=e[t], pen=QPen(Qt.red))

def plotLine(name, **data):
	for line in traceWidget.getViewBox().addedItems:
		if line.__name__ == name:
			line.setData(**data)
			return
	traceWidget.addItem(pg.PlotDataItem(**data), name=name, export=False)

def getPuffIndices(x, y, frame=-1):
	'''
	returns puff ID, Puff Group ID as a tuple, or just the group id if frame is not given
	'''
	if frame == -1:
		for i, group in enumerate(K['puffGroups']):
			if np.linalg.norm((round(group.x) - x, round(group.y) - y)) < .7:
				return -1, i
	else:
		for idx in K['puffGroups'][currentGroup].puff_ids:
			if K['puffs'][idx].t_start <= frame <= K['puffs'][idx].t_end:
				return idx, currentGroup
		idx = (frame, y, x)
	return -1, -1

def plotPuff(puff):
	puffTraceWidget.clear()
	puffTraceWidget.addItem(pg.PlotDataItem(x=puff.x_full, y=puff.f_full, pen=pg.mkPen('w')), name='full')
	puffTraceWidget.addItem(pg.PlotDataItem(x=puff.x_full, y=puff.weak_full, pen=pg.mkPen('r')), name='weak')
	puffTraceWidget.addItem(pg.PlotDataItem(x=puff.x_full, y=puff.strong_full, pen=pg.mkPen('g')), name='strong')

	puffTraceWidget.addItem(pg.InfiniteLine(pos=puff.baseline, angle=0, pen = .3), export=False)
	puffTraceWidget.addItem(pg.InfiniteLine(pos=puff.thresh20, angle=0, pen = .3), export=False)
	puffTraceWidget.addItem(pg.InfiniteLine(pos=puff.thresh50, angle=0, pen = .3), export=False)
	puffTraceWidget.addItem(pg.InfiniteLine(pos=puff.thresh80, angle=0, pen = .3), export=False)
	puffTraceWidget.addItem(pg.InfiniteLine(pos=puff.baseline+puff.amplitude, angle=0, pen =.3), export=False)
	puffTraceWidget.addItem(pg.ScatterPlotItem(pos = [(puff.t_start, puff.baseline)], symbol='o', symbolSize=10, pen=1), export=False)
	puffTraceWidget.addItem(pg.ScatterPlotItem(pos = [(puff.t_end, puff.weak_full[puff.t_end])], symbol='o', symbolSize=10, pen=1), export=False)
	puffTraceWidget.addItem(pg.ScatterPlotItem(pos = [(puff.t_peak, puff.baseline+puff.amplitude)], symbol='o', symbolSize=10, pen=1), export=False)

	puffTraceWidget.setXRange(puff.xdata[0], puff.xdata[-1])

def plot_puffs_on_image(): # points are plotted as (y, x)
	ps = map(lambda g: (g.x, g.y), K['puffGroups'])
	x, y = np.transpose(ps)
	data = {'x': x, 'y': y, 'size': 2, 'pen': QPen(Qt.yellow), \
		'name': 'Puff Group', 'pxMode': False, 'symbol': 's'}
	if not hasattr(imageWidget, 'puffScatter'):
		imageWidget.puffScatter = imageWidget.addItem(pg.ScatterPlotItem(**data), name='puffScatter')
	else:
		imageWidget.puffScatter.setData(**data)

def update_view(d):
	global K
	if d['viewDropDown'] == 'Average Frame':
		imageWidget.setImage(np.average(K['image'], 0), signal=False)
	elif d['viewDropDown'] == 'Video':
		imageWidget.setImage(K['image'], signal=False)
	K['roi_effect'] = d['ROI Effect']
	#update with group radius
	K['puffGroups'] = groupPuffs(K['puffs'], d['groupRadius'])
	plot_puffs_on_image()


def get_nearest_puffgroup(idx=0):
	x, y = imageWidget.location.getPosition() #image is transposed
	order = sorted([i for i in range(len(K['puffGroups']))],\
	 key=lambda i: np.linalg.norm((x - K['puffGroups'][i].x, y - K['puffGroups'][i].y)))
	i = order[idx]
	return i

def circleGroup(ids=None):
	global currentGroup, K
	print(ids)
	groupsToMerge = sorted(ids)
	if len(groupsToMerge) > 1:
		puff_ids=[]
		for i in groupsToMerge[::-1]:
			puff_ids.extend(K['puffGroups'][i].puff_ids)
			del K['puffGroups'][i]
		(x,y)=getpuffGroupOrigin(K['puffs'], puff_ids)
		K['puffGroups'].append(PuffGroup(puff_ids,x,y))
		currentGroup=len(K['puffGroups'])-1
		plot_puffs_on_image()
		plot_traces()

def onKeyPress(ev):
	global currentGroup
	key = ev.key()
	x, y = imageWidget.location.getPosition() # image is transposed
	if key==16777220: #enter
		if len(K['puffGroups']) == 1:
			print('Last Puff')
		else:
			circleGroup(ids=(currentGroup, get_nearest_puffgroup(idx=1)))
	elif key in (16777235, 16777236):
		if K['puffGroups']:
			if currentGroup >= 0:
				currentGroup+=1
				if currentGroup>=len(K['puffGroups']):
					currentGroup=0
			else:
				currentGroup = get_nearest_puffgroup()
	elif key in (16777234, 16777237):
		if K['puffGroups']:
			if currentGroup >= 0:
				currentGroup-=1
				if currentGroup<0:
					currentGroup=len(K['puffGroups']) - 1
			else:
				currentGroup = get_nearest_puffgroup()
	elif key in (68, 16777223):
		if K['puffGroups']:
			del K['puffGroups'][currentGroup]
			if currentGroup>=len(K['puffGroups']):
				currentGroup=0
			plot_puffs_on_image()
	imageWidget.location.setPos(QPoint(round(K['puffGroups'][currentGroup].x), round(K['puffGroups'][currentGroup].y)))

def trace_context(event):
	global currentGroup
	pt = event.screenPos()
	frame = round(traceWidget.getViewBox().mapToView(event.pos()).x())
	i = 0
	y, x = imageWidget.location.getPosition()
	puff_id, group_id = getPuffIndices(x, y, frame=frame)
	if group_id >= 0:
		currentGroup = group_id
	else:
		traceWidget.getViewBox().menu.exec_(pt)
		return
	for i in range(len(traceWidget.getViewBox().addedItems))[::-1]:
		it = traceWidget.getViewBox().addedItems[i]
		if it.__name__.startswith('B') and frame in it.getData()[0]:
			delete_puff(currentGroup, puff_id)
			return
		i += 1

def imageClick(ev):
	if not imageWidget.view.creatingROI:
		imageWidget.location.setPos(QPoint(imageWidget.view.mouse_x, imageWidget.view.mouse_y))
	ev.accept()
	ROIViewBox.mousePressEvent(imageWidget.view, ev)


def imageDrag(ev):
	if imageWidget.view.creatingROI:
		ROIViewBox.mouseDragEvent(imageWidget.view, ev)
	else:
		ev.accept()
		imageWidget.location.setPos(QPoint(imageWidget.view.mouse_x, imageWidget.view.mouse_y))

def traceClick(ev):
	if ev.button() == Qt.RightButton:
		trace_context(ev)
		ev.accept()
	else:
		ev.accept()
		pg.ViewBox.mousePressEvent(traceWidget.getViewBox(), ev)
		traceWidget.mouseLine.setValue(traceWidget.getViewBox().mapToView(ev.pos()).x())

def traceDrag(ev):
	ev.accept()
	traceWidget.mouseLine.setValue(traceWidget.getViewBox().mapToView(ev.pos()).x())

def mouseLine_moved(i_line):
	x = round(traceWidget.mouseLine.value())
	global currentGroup
	#traceWidget.mouseLine.setValue(x)
	pY, pX = imageWidget.location.getPosition()
	pY, pX = int(pY), int(pX)
	puff_id, currentGroup = getPuffIndices(pX, pY, frame=x)
	p = None
	if puff_id >= 0:
		ttl='Puff: %s' % (puff_id+1)
		if currentGroup >= 0 and currentGroup < len(K['puffGroups']):
			p=K['puffs'][puff_id]
			ttl='%s  Group: %i, Origin: (%f,%f)' % (ttl, currentGroup+1, p.y, p.x)
			plotPuff(p)
		else:
			ttl='%s (deleted)' % ttl
		puffTraceWidget.setTitle(ttl)

def roiCreated(roi):
	global K
	ids = [i for i, p in enumerate(K['puffGroups']) if roi.contains(QPoint(p.x, p.y))]
	if K['roi_effect'] == 'Delete Puffs':
		for i in sorted(ids)[::-1]:
			del K['puffGroups'][i]
	elif K['roi_effect'] == 'Join Puffs':
		circleGroup(ids)
	imageWidget.view.removeItem(roi)

def load_plots():
	global traceWidget, imageWidget, K
	traceWidget.scene().sigMouseMoved.disconnect(traceMouseMoved)

	traceWidget.removeItem(averageLine)
	K['puffGroups'] = groupPuffs(K['puffs'])
	mt, my, mx = np.shape(K['image'])
	imageWidget.location = PixelSelector(max_bounds = [my+1, mx+1]) # Image is transposed for user
	imageWidget.addItem(imageWidget.location, name='Mouse Position')
	imageWidget.view.mousePressEvent = imageClick
	imageWidget.view.mouseDragEvent = imageDrag
	imageWidget.location.removeHandle(0)
	imageWidget.location.sigRegionChanged.connect(plot_traces)
	imageWidget.view.roiCreated.connect(roiCreated)

	K['roi_effect'] = 'Delete Puffs'
	buttonWidget = OptionsWidget('', [{'key': 'viewDropDown', 'name':'Image Background', \
		'value': ['Video', 'Average Frame']}, {'key': 'ROI Effect', 'value': ['Delete Puffs', 'Join Puffs']}, \
		{'key': 'button', 'name': 'Restart', 'action': reset},\
		{'key': 'groupRadius', 'name': 'Group Radius', 'value': 1}])
	buttonWidget.valueChanged.connect(update_view)
	win.addWidget(buttonWidget, name='Buttons', where=('bottom', imageDock), size=(9, 1))

	traceWidget.mouseLine = pg.InfiniteLine(0, pen=QPen(Qt.blue), bounds = (0, len(K['image'])))
	traceWidget.addItem(traceWidget.mouseLine, name='Frame Position', export=False)
	traceWidget.mouseLine.sigPositionChanged.connect(mouseLine_moved)
	traceWidget.getViewBox().mousePressEvent = traceClick
	traceWidget.getViewBox().mouseDragEvent = traceDrag

	traceWidget.keyPressEvent = onKeyPress
	imageWidget.keyPressEvent = onKeyPress
	puffTraceWidget.keyPressEvent = onKeyPress

	plot_puffs_on_image()
	plot_traces()

def delete_puff(group_id, puff_id):
	global currentGroup, K
	new_ids = [i for i in K['puffGroups'][group_id].puff_ids if i != puff_id]
	if not new_ids:
		del K['puffGroups'][group_id]
		currentGroup = group_id-1
		if currentGroup < 0:
			currentGroup = len(K['puffGroups']) - 1
		imageWidget.location.setPos(QPoint(round(K['puffGroups'][currentGroup].x), round(K['puffGroups'][currentGroup].y)))
	else:
		[ynew,xnew]=getpuffGroupOrigin(K['puffs'], new_ids)
		K['puffGroups'][group_id] = PuffGroup(new_ids, xnew, ynew)
		imageWidget.location.setPos(QPoint(xnew, ynew))
	plot_puffs_on_image()

def reset():
	traceWidget.getViewBox().clear()
	imageWidget.view.clear()
	traceWidget.scene().sigMouseMoved.connect(traceMouseMoved)
	traceWidget.addItem(averageLine, name='Average Trace')
	imageWidget.open_image()
	win.moveDock(opsDock, 'above', puffDock)

def start_flika():
	global K, currentThread
	currentThread = FLIKA_Thread(K)
	currentThread.done.connect(load_plots)
	currentThread.step.connect(imageWidget.setImage)
	beginAction.setEnabled(False)
	currentThread.start()

class FLIKA_Thread(QThread):
	step = Signal(np.ndarray)
	done = Signal(dict)
	def __init__(self, K):
		super(FLIKA_Thread, self).__init__()
		self.K = K

	def run(self):
		spatial = spatially_filter(self.K['image'], self.K['sigmaXY'])
		self.step.emit(spatial)
		temporal = temporally_filter(spatial, self.K['tfiltN'], K['flash_frames'])
		self.step.emit(temporal)
		self.K['ratioed'], self.K['ratioed2'] = calculate_df_f(spatial, temporal)
		m_filt_norm = calculate_threshold(spatial, temporal)
		self.step.emit(m_filt_norm)
		del spatial
		self.K['puffBoolean'] = detect_puffing_pixels(m_filt_norm, self.K)
		#self.step.emit(self.K['puffBoolean'] * m_filt_norm)
		del m_filt_norm
		self.K['puff_idx'] = group_into_puffs(K['puffBoolean'], self.K)
		self.K['puffs'] = analyze_puffs(self.K['puff_idx'], self.K, temporal)
		self.step.emit(self.K['image'])
		self.done.emit(self.K)

def imageMouseMoved(pos):
	pos = imageWidget.getImageItem().mapFromScene(pos)
	s = 'Mouse at (%s, %s)' % (int(pos.x()), int(pos.y()))
	if imageWidget.loaded():
		t, y, x = K['image'].shape
		if QRectF(0, 0, y - 1, x - 1).contains(pos):
			frame = K['image'][imageWidget.currentIndex]
			val = frame[int(pos.x()), int(pos.y())]
			s += ' - Value: %.3f - Frame Min/Max = %.2f/%.2f' % (val, np.min(frame), np.max(frame))
		else:
			s += ' - Out of bounds'
	win.statusBar().showMessage(s)

def traceMouseMoved(pos):
	pos = traceWidget.getPlotItem().mapToView(pos)
	pos = QPoint(int(pos.x()), int(pos.y()))
	s = 'Mouse at Frame %d' % pos.x()
	if imageWidget.loaded():
		t = averageLine.getData()[1]
		if pos.x() >= 0 and pos.x() < len(t):
			s += ' - Value: %.3f - Trace Min/Max = %.2f/%.2f' % (t[pos.x()], np.min(t), np.max(t))
			imageWidget.setCurrentIndex(pos.x())
		else:
			s += ' - Outside of Frame Range'
	else:
		s += ' - No image loaded.'
	win.statusBar().showMessage(s)

def saveData():
		dirname = str(QFileDialog.getExistingDirectory(win, "Select Directory to save all info to..."))

		cons = ['crop_frames', 'flashExists', 'black_level_val', 'puffDetectionMethod', 'sigmaXY', 'tfiltN', 'threshold_constant', \
			'min_rise_time', 'max_rise_time', 'dilateXY', 'dilateT', 'minPuffSize', 'padding', 'wfiltN', 'sfiltN', 'max_fall_time', 'RotatedFit', 'maxSigmaForGaussianFit']
		with open(os.path.join(dirname, 'Contants.txt'), 'w') as outf:
			outf.write('\t'.join(cons))
			for i, con in enumerate(cons):
				if type(K[con]) not in (str, bool, tuple, list):
					outf.write(float(K[con]))
				else:
					outf.write(str(K[con]))
				outf.write('\t')

		with open(os.path.join(dirname, 'Puff Data.txt', 'w')) as outf:
			pass
		# Save all the information about each group and puff
		data = []
		descs=['Group #','Group x','Group y','N Events','Max Amp','x','y','t_peak','Amplitude','sigmax','sigmay','angle','r20','r50','r80','r100','f80','f50','f20','f0']
		if puffs:
			for i in range(len(puffGroups)):
				puffGroup = puffGroups[i]
				for idx in puffGroup.puff_ids:
					p = puffs[idx]
					data.append(dict(zip(descs, [i+1, puffGroup.x, puffGroup.y, "", "", p.x, p.y, p.t_peak, p.amplitude, p.sigmax, p.sigmay, p.angle, p.r20, p.r50, p.r80, p.r100, p.f80, p.f50, p.f20, p.f0])))
			freq = 1
			maxAmp = data[0]['Amplitude']
			groupN = data[0]['Group #']
			data = [d for d in data if d[0] != np.nan]
			for i in range(1, len(data)):
				if data[i]['Group #'] == groupN:
					freq += 1
					maxAmp = max(maxAmp, data[i][8])
					data[i]['Group x'] = ""
					data[i]['Group y'] = ""
				else:
					data[i-freq]['N Events'] = freq
					data[i-freq]['Max Amp'] = maxAmp
					groupN = data[i]['Group #']
					freq = 1
					maxAmp = data[i]['Amplitude']
			num =len(data)
			data[num-freq]['N Events'] = freq
			data[num-freq]['Max Amp'] = maxAmp


			export_dict_list(data, os.path.join(dirname, 'Puff Data.txt'), order=descs)

			# now save the group traces
			data=np.zeros((np.size(K["Ratioed"],0), len(puffGroups)))
			for i in range(len(puffGroups)):
				puffGroup=puffGroups[i]
				x=round(puffGroup.x)
				y=round(puffGroup.y)
				data[:, i]=np.squeeze(K['Ratioed'][:, y, x])

			with open(os.path.join(dirname, 'Group Traces.txt'), 'w') as outf:
				out.write("\t".join(["Group %d" % i for i in range(np.size(data, 1))]))
				for i in range(np.size(data, 0)):
					outf.write("\t".join([(str(data[i][gnum]) if not np.isnan(data[i][gnum]) else "") for gnum in range(np.size(data, 1))]))

imageWidget.imageChanged.connect(imageChanged)
traceWidget.scene().sigMouseMoved.connect(traceMouseMoved)
imageWidget.scene.sigMouseMoved.connect(imageMouseMoved)
traceWidget.load_file = imageWidget.open_image

fileMenu = win.menuBar().addMenu('&File')
fileMenu.addAction(QAction('&Open Image', fileMenu, triggered = lambda : imageWidget.open_image()))
fileMenu.addAction(QAction('&Close', fileMenu, triggered = win.close))
fileMenu.addAction(QAction('&Reset', fileMenu, triggered = reset))
analyzeMenu = win.menuBar().addMenu("&FLIKA")
beginAction = QAction('&Analyze', analyzeMenu, triggered = start_flika)
analyzeMenu.addAction(beginAction)
exportMenu = win.menuBar().addMenu('&Export')
exportMenu.addAction(QAction('&Save Data', exportMenu, triggered=saveData))
beginAction.setEnabled(False)
traceDock = win.addWidget(traceWidget, name="Full Trace", size=(5, 4))
puffDock = win.addWidget(puffTraceWidget, name='Puff Trace', size=(3, 1), where=('bottom', traceDock))
make_opts()
opsDock = win.addWidget(opsWidget, name='Analysis Options', size=(2, 1), where=('above', puffDock))
imageDock = win.addWidget(imageWidget, name='Tiff File', size=(2, 1), where=('right', puffDock))

win.show()
app.exec_()
