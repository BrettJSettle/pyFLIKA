"""
@author: Brett Settle
@Department: UCI Neurobiology and Behavioral Science
@Lab: Parker Lab
@Date: August 6, 2015
"""
from math import pi, cos, sin
from scipy.optimize import curve_fit, basinhopping
import numpy as np
from leastsqbound import leastsqbound

def sind(num):
    return sin(num * pi/180.)
def cosd(num):
    return cos(num * pi/180.)

def fitGaussian(I, r, maxSigmaForGaussianFit, rotatedFit):
    # takes an nxm matrix and returns an nxm matrix which is the gaussian fit
    # of the first.  c is a list of parameters [zoffset, scale, xorigin, yorigin, sigma]
    I[I<0]=0
    startx=round(np.mean(r[0:2]))
    starty=round(np.mean(r[2:4]))
    center=I[r[0]:r[1],r[2]:r[3]]
    scale=np.max(center.ravel())
    if scale>0:
        I=np.divide(I,scale)
    n,m=np.shape(I)
    Y,X=np.meshgrid(range(1, m+1),range(1, n+1))
    x = np.zeros((np.size(X), 2))
    x[:, 0] = X.ravel('F')
    x[:, 1] = Y.ravel('F')
    I_col=I.ravel('F')
    if not rotatedFit: # fun = e ^ (-())
        def fun(c, x, minus=True):
            if minus==False:
                return np.exp(-((x[:, 0]-c[0])/c[2])**2-((x[:, 1]-c[1])/c[2])**2)
            return I_col - np.exp(-((x[:, 0]-c[0])/c[2])**2-((x[:, 1]-c[1])/c[2])**2)
    #%     [  xorigin,  yorigin, sigma]
        c0=[startx,starty,m] #starting parameters
        lb=[r[0],r[2],1 ] #lower bounds on the parameters
        ub=[r[1],r[3],maxSigmaForGaussianFit] #upper bounds on the parameters
        for i in range(len(c0)):
            c0[i] = max(lb[i], c0[i])
            c0[i] = min(ub[i], c0[i])
        c = leastsqbound(fun, c0, args=(x,), bounds = zip(lb, ub), ftol = 1.0e-4, maxfev=100)

    else: #if we are fitting to a 2d rotating gaussian fit
        def fun(c, x, minus=True):
            if not minus:
                return np.exp(-(((x[:, 0] - c[0])*cosd(c[4])+(x[:, 1] - c[1])*sind(c[4]))/c[2])**2 - (((x[:, 0] - c[0])*sind(c[4])+(x[:, 1] - c[1])*cosd(c[4]))/c[3])**2)
            return I_col - np.exp(-(((x[:, 0] - c[0])*cosd(c[4])+(x[:, 1] - c[1])*sind(c[4]))/c[2])**2 - (((x[:, 0] - c[0])*sind(c[4])+(x[:, 1] - c[1])*cosd(c[4]))/c[3])**2)
    #%     [xorigin yorigin sigmax sigmay angle]
        c0=[startx, starty, m, m, 90] #starting parameters
        lb=[r[0], r[2], 1, 1, 0] #lower bounds on the parameters
        ub=[r[1], r[3], maxSigmaForGaussianFit, maxSigmaForGaussianFit, 180] #upper bounds on the parameters
        for i in range(len(c0)):
            c0[i] = max(lb[i], c0[i])
            c0[i] = min(ub[i], c0[i])
        c = leastsqbound(fun, c0, args=(x,), bounds = zip(lb, ub), ftol = 1.0e-6, maxfev=100)
    c = c[0]
    Ifit=fun(c, np.double(x), minus=False)
    Ifit=np.reshape(Ifit,(n, m), 'F')
    graph_args = [I, Ifit, X, Y]

    return c, graph_args
