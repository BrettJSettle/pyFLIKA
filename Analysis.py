"""
@author: Brett Settle
@Department: UCI Neurobiology and Behavioral Science
@Lab: Parker Lab
@Date: August 6, 2015
"""
from scipy.signal import convolve, filtfilt, butter
import numpy as np
from collections import namedtuple
from fitGaussian import fitGaussian
from scipy.ndimage.measurements import label
from scipy.ndimage.morphology import grey_dilation
from scipy.ndimage.filters import gaussian_filter
from math import sin, cos, pi, sqrt, ceil, floor

PuffGroup = namedtuple("PuffGroup", "puff_ids x y")

def locate_flash(arr):
    arr = np.array(arr)
    if arr.ndim > 1:
        arr = np.median(arr, (2, 1))
    medianSlopes = np.diff(arr, 1, 0)
    medianSlopes = np.insert(medianSlopes, 0, 0)
    median_std = np.std(medianSlopes)
    start=np.where(medianSlopes > 2 * median_std)[0]
    start=[frame for frame in start if frame > 10] #needs to be after the 10th frame
    try:
        start=start[0]
    except:
        return 0, 0

    end=np.where(medianSlopes < -2 * median_std)[0]
    end=[frame for frame in end if frame > start]
    try:
        end=end[-1]
    except:
        return 0, 0

    del medianSlopes, median_std

    start=start-5
    end=end+6
    return start, end

def spatially_filter(arr, sigmaXY):
    print('Spatially filtering')
    spatial = np.zeros_like(arr)
    if sigmaXY > 0:
        gaussianKernalSize=int(np.ceil((4*sigmaXY+1)/2)*2+1)
        for t in np.arange(np.size(arr, 0)):
            #if len(flash) == 0 or not (flash[0] <= t <= flash[1]):
            #spatial[t]=cv2.GaussianBlur(arr[t],(gaussianKernalSize,gaussianKernalSize),sigmaXY)
            spatial[t] = gaussian_filter(arr[t], sigmaXY)
    else:
        spatial = np.copy(arr)
    print('Spatially Filtered')
    return spatial

def temporally_filter(spatial, tfiltN, flash):
    print("Temporally Filtering")
    padlen=0
    if tfiltN[2]==1: #%if only high pass temporal filter
        [b,a]= butter(tfiltN[0],tfiltN[1],btype='highpass')
        padlen=3
    else: #%if band pass temporal filter
        [b,a]=butter(tfiltN[0],[tfiltN[1],tfiltN[2]], btype='bandpass')
        padlen=6
    temporal=np.array(spatial, dtype = np.float32)

    if len(flash) == 2:
        temporal[flash[0]:flash[1]+1,:,:]=np.mean(spatial[4:flash[0]+1, :, :], 0)
        bg=np.mean(temporal[:flash[0] + 1, :, :],0)
    else:
        bg = np.min(temporal, 0)

    for t in np.arange(len(temporal)):
        temporal[t, :, :] -= bg

    #temporal[0, :, :] = 0
    p = 0.0
    for y in range(np.size(temporal, 1)):
        for x in range(np.size(temporal, 2)):
            temporal[:, y, x]=filtfilt(b,a, temporal[:, y, x], padlen=padlen)
        if (100 * y / np.size(temporal, 1)) > p:
            p = (100 * y / np.size(temporal, 1))
            print("  %d%%\r" % (100 * y / np.size(temporal, 1))),
    print('Temporally Filtered')
    return temporal

def calculate_df_f(spatial, temporal):
    # calculate deltaf/f  by  Dividing Each Pixel by Baseline %%%%%%
    print("Calculating dF / F0")
    #if len(flash) == 2:
    #    normalize_M=np.median(spatial[4:flash[0]+1, :, :],0)
    #else:
    normalize_M=np.median(spatial,0) * .75 # 3/4 of full median as normal, so it's not too high
    ratioed=np.ones_like(spatial)
    ratioed2=np.zeros_like(temporal)
    np.seterr(divide='ignore', invalid='ignore')

    for t in np.arange(np.size(spatial,0)):
        #if len(flash) == 2 and flash[0] <= t <= flash[1]:
        #    continue
        ratioed[t] = spatial[t] / normalize_M
        ratioed2[t] = spatial[t] / normalize_M

    print("dF/F0 calculated")
    return ratioed, ratioed2


def calculate_threshold(spatial, temporal):
    #% Calculate Threshold for each Pixel
    print("Calculating Threshold")
    [b,a]=butter(1,.2,'highpass')
    #if len(flash) == 2:
    #    B_filt=np.array(spatial[4:flash[0]-1, :, :], dtype=np.float32)
    #else:
    B_filt=np.array(spatial, dtype=np.float32)
    mz, my, mx = np.shape(spatial)
    B_filt -= B_filt[0]
    p = 0.0
    for y in np.arange(my):
        for x in np.arange(mx):
            B_filt[:, y, x] = filtfilt(b, a, B_filt[:, y, x], padlen=3)
        if (100 * y / my) > p:
            p = (100 * y / my)
            print ("  %d%%\r" % p),

    B_std=np.std(B_filt, 0, ddof=1)

    m_filt_norm=np.zeros_like(temporal, dtype=np.float32)
    for t in np.arange(np.size(temporal,0)):
        m_filt_norm[t, :, :] = temporal[t, :, :] / B_std
    print('Threshold Calculated')
    return m_filt_norm

def detect_puffing_pixels(m_filt_norm, opts):
    #Subtract pixel intensities from different frames and compare to threshold
    mz, my, mx = np.shape(m_filt_norm)

    print("Detecting Puffing Pixels")
    f_pad=5 #if you change this constant, change the corresponding constant in analyzePuff
    if opts['puffDetectionMethod'] == 'slope':
        puffBoolean=np.zeros((mz, my, mx), dtype = bool) # return this
        puff_idx=[]
        for dt in np.unique(np.round(np.linspace(opts['min_rise_time'], opts['max_rise_time'], 10, endpoint=True))): #min_rise_time:max_rise_time
            # try 10 shifts between min and max to find puffing pixels
            other = np.subtract(m_filt_norm[dt:, :, :].ravel(), m_filt_norm[:-dt, :, :].ravel())
            puff_idx.extend([i + dt*my*mx for i in np.where(other > opts['threshold_constant'])[0]])
            print("  %d%%\r" % (100. * dt / opts['max_rise_time'])),
        print("Puffing Pixels Detected")
        puff_idx=np.double(puff_idx)
        unqP=np.unique(puff_idx)
        countElP= np.histogram(puff_idx, unqP)[0]
        puff_idx=np.array(unqP[countElP>=2], dtype=np.uint32)
        puff_idx = unravel_indices(puff_idx, np.shape(puffBoolean))
        print("Puffing Pixels Located")
        for t, y, x in puff_idx:
            puffBoolean[t, y, x]=1
        puffBoolean[:opts['max_rise_time'] + 1, :, :]=0
        puffBoolean[-f_pad:, :, :]=0
    elif opts['puffDetectionMethod'] == 'threshold':
        puffBoolean=opts['m_filt_norm']>opts['threshold_constant']
        puffBoolean[:opts['max_rise_time'] + 1, :, :]=0
        puffBoolean[-opts['max_fall_time']:, :, :]=0
        puffBoolean[opts['flash_frames'][0]:opts['flash_frames'][1] + 1, :, :]=0
        puffBoolean[-f_pad:, :, :]=0
    return np.array(puffBoolean)

def group_into_puffs(puffBoolean, opts):
    '''
    takes the puffBoolean image and the options dictionary, returns the puff ids and their locations
    '''
    print('Grouping Into Puffs')
    puffBooleanDilate = grey_dilation(puffBoolean, (opts['dilateT'], opts['dilateXY'], opts['dilateXY']))
    results=bwconncomp(puffBooleanDilate, puffBoolean) #join all puffs together
    del puffBooleanDilate
    print("%d puffs detected" % results.NumObjects)
    puff_idx = [puff for puff in results.PixelIdxList if len(puff) >= opts['minPuffSize']]
    print('%d puffs remaining because some were too small.' % np.size(puff_idx, 0))
    return puff_idx

def analyze_puffs(puff_idx, K, temporal):
    print("Analyzing puffs")
    pcompleted=-1
    puffs=[]
    puff_count = np.size(puff_idx, 0)
    for j in np.arange(puff_count):
        puffs.append(Puff(puff_idx[j], temporal))
        puffs[-1].get_bounds(K['padding'])
        puffs[-1].analyze(K)
        if pcompleted<floor(100 * (j/puff_count)):
            pcompleted=floor(100 * (j/puff_count))
            print("  %d%%\r" % pcompleted),

    print("Puffs Analyzed")
    return puffs


def bwconncomp(dilated, boolean):
    Dilation = namedtuple('Dilation', 'Connectivity ImageSize NumObjects PixelIdxList')
    size = dilated.shape
    labeled_frame, objects = label(dilated, np.ones((3, 3, 3)))
    labeled_frame *= boolean
    pixel_list = []
    for i in np.arange(1, objects + 1):
        raveled = labeled_frame.ravel()
        pos = np.where(raveled==i)[0]
        puff = unravel_indices(pos, np.shape(labeled_frame))
        pixel_list.append(puff)
        print("On Puff %d of %d          \r" % (i, objects)),

    del labeled_frame

    return Dilation(26 if len(size) == 3 else 8, size, objects, pixel_list) # list of [t, y, x] lists

def unravel_indices(arr, shape):
    return np.array(zip(*np.unravel_index(arr, shape)))

def ravel_indices(arr, shape):
    return np.array([np.ravel_multi_index(arr[i], shape) for i in range(len(arr))])

def fit_to_polynomial(self, trace, baseline = 0):
    trace = np.copy(trace)
    x=np.arange(len(trace))
    np.warnings.simplefilter('ignore', np.RankWarning)
    poly=np.poly1d(np.polyfit(x,trace,20))
    info = {}
    info['x'] = x
    info['yfit'] = poly(x)
    info['amplitude']=np.max(ftrace)
    info['baseline'] = baseline

    info['Rise 20'] = -1
    info['Rise 50'] = -1
    info['Rise 80'] = -1
    info['Rise 100'] = -1
    info['Fall 80'] = -1
    info['Fall 50'] = -1
    info['Fall 20'] = -1
    try:
        info['Threshold 20']=(info['amplitude'] - info['baseline'])*.2 + info['baseline']
        info['Threshold 50']=(info['amplitude'] - info['baseline'])*.5 + info['baseline']
        info['Threshold 80']=(info['amplitude'] - info['baseline'])*.8 + info['baseline']
        pt = 'Rise 20'
        info['Rise 20']=np.argwhere(info['yfit']>info['Threshold 20'])[0][0]
        pt = 'Rise 50'
        info['Rise 50']=np.argwhere(info['yfit']>info['Threshold 50'])[0][0]
        pt = 'Rise 80'
        info['Rise 80']=np.argwhere(info['yfit']>info['Threshold 80'])[0][0]
        pt = 'Rise 100'
        info['Rise 100']=np.argmax(info['yfit'])
        pt = 'Fall 80'
        tmp=np.squeeze(np.argwhere(info['yfit']<info['Threshold 80']))
        info['Fall 80']=tmp[tmp>info['Rise 100']][0]
        pt = 'Fall 50'
        tmp=np.squeeze(np.argwhere(info['yfit']<info['Threshold 50']))
        info['Fall 50']=tmp[tmp>info['Fall 80']][0]
        pt = 'Fall 20'
        tmp=np.squeeze(np.argwhere(info['yfit']<info['Threshold 20']))
        info['Fall 20']=tmp[tmp>info['Fall 50']][0]
    except Exception as e:
        print("Analysis failed when calculating the %s point.")
    return info

def groupPuffs(puffs, radius=1):
    puffGroups = []
    radii=np.zeros((len(puffs),len(puffs)))
    radii.fill(200)
    for i in range(len(puffs)):
        for j in range(i, len(puffs)):
            radii[i,j]=sqrt((puffs[i].x-puffs[j].x)**2+(puffs[i].y-puffs[j].y)**2)
    val = max([sqrt(2),radius])
    near=radii<=val
    for i in range(len(puffs)-1, -1, -1): # for each puff, starting at the last
        for j in range(i-1, -1, -1):    # for each puffs distance from it boolean
            if any(near[i,:] & near[j,:]): # if booleans are both true
                near[j, :]=near[i,:] | near[j,:] # join the group
                near = np.delete(near, i, 0)
                break
    for i in range(np.size(near,0)):
        puff_ids=np.where(near[i,:])[0]
        ynew, xnew=getpuffGroupOrigin(puffs, puff_ids)
        puffGroups.append(PuffGroup(puff_ids, xnew, ynew))
    return puffGroups

def getpuffGroupOrigin(puffs, puff_ids):
    x=y=0
    kk = len(puff_ids)
    for k in puff_ids:
        x+=puffs[k].x
        y+=puffs[k].y
    x/=kk
    y/=kk
    return x, y

def findextrema(v):
    v_prime=np.concatenate(([0], np.diff(v)))
    zero_crossing=np.concatenate((np.diff(np.sign(v_prime)),[0]))
    maxima=[1 if v < -1 else 0 for v in zero_crossing]
    minima=[1 if v > 1 else 0 for v in zero_crossing]
    if np.sign(v_prime[-1])>0:
        maxima[-1]=1
    else:
        minima[-1]=1;
    if np.sign(v_prime[0])>0:
        maxima[0]=1
    else:
        minima[0]=1
    maxima=np.nonzero(maxima)[0]
    minima=np.nonzero(minima)[0]
    return maxima, minima

class Puff:
    def __init__(self, pixel_idx, temporal):
        self.puff_idx=pixel_idx
        self.temporal = temporal
        # get puff pixel coordinates and the video to filter

    def get_bounds(self, padding):

        self.p = [padding, padding, padding, padding]
        mt, my, mx = np.shape(self.temporal)
        ts, ys, xs = zip(*self.puff_idx)
        self.before=min(ts)
        self.after=max(ts)
        self.top=min(ys)-self.p[0]
        self.bottom=max(ys)+self.p[1]
        self.left=min(xs)-self.p[2]
        self.right=max(xs)+self.p[3]
        if self.top < 0:
            self.top=0
            self.p[0]=min(ys)
        if self.left < 0:
            self.left=0
            self.p[2]=min(xs)
        if self.bottom>my:
            self.bottom=my
            self.p[1]=my-max(ys)
        if self.right>mx:
            self.right=mx
            self.p[3]=mx-max(xs)

    def analyze(self, K):
        mt, my, mx = np.shape(self.temporal)
        ts, ys, xs = zip(*self.puff_idx)
        self.cropped = self.temporal[self.before:self.after+1,self.top:self.bottom+1,self.left:self.right+1]
        I=np.mean(self.temporal[self.before:self.after+1,self.top:self.bottom+1,self.left:self.right+1],0)
        # index is [Y, X]
        rng=[self.p[0], np.size(I,0)-self.p[1], self.p[2], np.size(I,1)-self.p[3]] # top bottom left right
        c, self.graph_args = fitGaussian(I,rng, K['maxSigmaForGaussianFit'], K['RotatedFit']) #causes math to fail
        origin=(c[0], c[1])
        self.y=origin[0]+self.top
        self.x=origin[1]+self.left

        if K['RotatedFit']:
            self.sigmay=c[2]
            self.sigmax=c[3]
            self.angle=c[4]
        else:
            self.sigmax=np.nan
            self.sigmay=np.nan
            self.angle=np.nan

        self.before = self.before - K['max_rise_time']
        self.after = self.after + K['max_fall_time']
        b1,a1=butter(K['sfiltN'][0],K['sfiltN'][1],btype='lowpass') #for fStrong
        b2,a2=butter(K['wfiltN'][0],K['wfiltN'][1],btype='lowpass') #for fWeak

        f_pad=5 #if you change this constant, change the corresponding constant in detect_puffing_pixels
        self.after=min(mt-f_pad, self.after)
        self.before=max(self.before, 1+f_pad)
        f=np.squeeze(K['ratioed'][self.before-f_pad:self.after+f_pad+1, round(self.y), round(self.x)])
        self.f_full = np.squeeze(K['ratioed'][:, round(self.y), round(self.x)])
        self.x_full = np.arange(mt)
        self.strong_full = filtfilt(b1, a1, self.f_full)
        self.weak_full = filtfilt(b2, a2, self.f_full)


        fStrong=filtfilt(b1,a1,f)
        fWeak = filtfilt(b2,a2,f)
        f=f[f_pad:-f_pad]
        fStrong=fStrong[f_pad:-f_pad]
        fWeak=fWeak[f_pad:-f_pad]
        self.xdata = np.arange(self.before, self.after+1)

        # use maximum slope to find the rising phase of the puff
        f_prime=np.concatenate(([0],np.diff(fStrong)))
        max_slope_idx = np.where(f_prime==np.max(f_prime[K['min_rise_time']:max(ts)-self.before]))[0][0]
        max_slope = f_prime[max_slope_idx]

        # use the local minimum immediately before the maximum slope as the start of the puff
        [maxima,minima]=findextrema(fStrong)
        try:
            t_start=max(minima[minima<max_slope_idx])
        except Exception:
            t_start=0
        baseline=fStrong[t_start]
        tmp=np.where(fWeak>baseline)[0]
        try:
            tmp = [val for val in tmp if t_start<=val<=max_slope_idx]
            t_start=tmp[0]
        except Exception:
            t_start=max_slope_idx

        # use the next local maxima of fStrong as the time of the peak
        t_peak=min(maxima[maxima>=max_slope_idx])

        del max_slope_idx, maxima

        # stop looking for a decrease in puff the next time the slope is a quarter the value of the maximum slope
        rising_times=np.where(f_prime>.25*max_slope)[0]

        try:
            t_end=min(rising_times[rising_times>t_peak])
            t_end=max(minima[minima<t_end])
        except Exception:
            t_end=len(f)-1

        #find a more accurate t_peak
        t_peak = np.where(fWeak==max(fWeak[t_start:t_end]))[0][0]
        self.f_peak=fWeak[t_peak] # get first value
        if baseline>self.f_peak:
            baseline=fWeak[t_start]

        self.baseline = baseline
        self.t_peak = t_peak
        self.t_end = t_end
        self.t_start = t_start
        self.f = f
        self.f_weak = fWeak
        self.f_strong = fStrong

        self.findRiseAndFall()

    def findRiseAndFall(self, end=True):
        self.f_peak = self.f_weak[self.t_peak]
        self.thresh20=self.baseline+(self.f_peak-self.baseline)*.2;
        self.thresh50=self.baseline+(self.f_peak-self.baseline)*.5;
        self.thresh80=self.baseline+(self.f_peak-self.baseline)*.8;
        tmp=np.where(self.f_weak>self.thresh20)[0]
        tmp=[val for val in tmp if val>=self.t_start and val<=self.t_peak]
        try:
            self.r20=tmp[0]
        except:
            print(tmp, self.t_start, self.t_peak)
        tmp=np.where(self.f_weak>self.thresh50)[0]
        tmp=[val for val in tmp if val>=self.t_start and val<=self.t_peak]
        try:
            self.r50=tmp[0]
        except:
            print(tmp, self.t_start, self.t_peak)
        tmp=np.where(self.f_weak>self.thresh80)[0]
        tmp=[val for val in tmp if val>=self.t_start and val<=self.t_peak]
        try:
            self.r80=tmp[0]
        except:
            print(tmp, self.t_start, self.t_peak)
        tmp=np.where(self.f_weak<self.thresh80)[0]
        tmp=[val for val in tmp if val>=self.t_peak]
        if not tmp:
            self.f80=np.nan
        else:
            self.f80=tmp[0]
        tmp=np.where(self.f_weak<self.thresh50)[0]
        tmp=[val for val in tmp if val>=self.t_peak]
        if not tmp:
            self.f50=np.nan
        else:
            self.f50=tmp[0]
        tmp=np.where(self.f_weak<self.thresh20)[0]
        tmp=[val for val in tmp if val>=self.t_peak]
        if not tmp:
            self.f20=np.nan
        else:
            self.f20=tmp[0]
        if end:
            tmp=np.where(self.f_weak<self.baseline)[0]
            tmp=[val for val in tmp if val>=self.t_peak]
            if not tmp:
                self.f0=np.nan
            else:
                self.f0=tmp[0]
            tmp=np.where(self.f_weak<self.baseline)[0]
            tmp=[val for val in tmp if val>=self.t_peak]
            if tmp and tmp[0]<self.t_end:
                self.t_end=tmp[0]
        else:
            self.f0 = self.t_end

        self.amplitude=self.f_peak-self.baseline
        self.r20-=self.t_start
        self.r50-=self.t_start
        self.r80-=self.t_start
        self.r100=self.t_peak-self.t_start
        self.f80-=self.t_peak
        self.f50-=self.t_peak
        self.f20-=self.t_peak
        self.f0-=self.t_peak
        self.t_peak+=self.before
        self.t_start+=self.before
        self.t_end+=self.before

    def set_start(self, x):
        self.t_start = max(0, min(x, self.t_peak - 1)) - self.before
        self.t_end -= self.before
        self.t_peak -= self.before
        self.findRiseAndFall(end=False)
        self.baseline = self.weak_full[self.t_start]
        self.set_peak(self.t_peak)

    def set_peak(self, x):
        if self.weak_full[x] <= self.weak_full[self.t_start] or self.weak_full[x] <= self.weak_full[self.t_end]:
            return
        self.t_peak = max(min(self.t_end-1, x), self.t_start+1)
        self.t_peak -= self.before
        self.t_end -= self.before
        self.t_start -= self.before
        self.findRiseAndFall(end=False)
        self.amplitude = self.weak_full[self.t_peak] - self.baseline
        self.analyze(False)

    def set_end(self, x):
        self.t_end = max(x, self.t_peak + 1)
        self.t_end -= self.before
        self.t_peak -= self.before
        self.t_start -= self.before
        self.findRiseAndFall(end=False)
        self.analyze(False)

    def __contains__(self, i):
        return i in self.puff_idx
