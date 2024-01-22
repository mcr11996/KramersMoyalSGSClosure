# -*- coding: utf-8 -*-
"""
Created on Sun May 7 01:18:16 2023

@author: Molly Ross
"""
import numpy as np


def transition_matrix(ts_dig,tau):
    """
    Finds the transition matrix of a digitized time series.

    Parameters
    ----------
    ts_dig : TYPE: numpy.ndarray with dtype=int64
        Time series digitized to the bin number with N bins.
    tau : TYPE: integer
        Number of time shifts to the next value.

    Returns
    -------
    TYPE: NxN numpy.ndarray with dtype=float
        Transition matrix describing probability of going from bin i (rows)
        to bin j (columns) after tau time steps.

    """
    n = 1+ max(ts_dig) # Number of states

    M = [[0]*n for _ in range(n)] # Initiate transition probability matrix

    for (i,j) in zip(ts_dig,ts_dig[tau:]):
        M[i][j] += 1
    
    # Convert to Fractions:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return np.asarray(M)

def findQ(ts_raw, lam_bins = 10, numPeriods=3,tau=20):
    """
    Calculate the Chi Square Statistic for transition probability of various
    sections of the data
    Ref: https://www.files.ethz.ch/isn/124233/kap1086.pdf (Section 1 p.8)

    Parameters
    ----------
    ts_raw : TYPE, numpy.ndarray with dtype=float.
        Raw time series to obtain the chi-sqaure statistic.
    num_bins : TYPE, integer.
        Number of bins to devide the data into.The default is 10.
    numPeriods : TYPE, integer
        Number of sections the time series will be divided into to calculate
        the Chi Square (Q) statistic. The default is 3.
    tau : TYPE, integer
        Time shift for calculating transition probability. The default is 20.

    Returns
    -------
    Q : TYPE, float.
        Chi-square statistic for the transition probabilities of the divided
        time series.

    """
    bins = np.linspace(np.min(ts_raw),np.max(ts_raw))
    ts_dig = np.digitize(ts_raw,bins)
    transMatrix = transition_matrix(ts_dig,tau = tau)
    splitMatrices = [] # Divide the time series into sub-series
    for i in range(numPeriods):
        splitMatrices.append(transition_matrix(ts_dig[int(i*len(ts_dig)/numPeriods):int((i+1)*len(ts_dig)/numPeriods)],tau=tau))

    Q=0
    for i in range(numPeriods):
        for j in range(lam_bins):
            rowsumt = sum(splitMatrices[i][j])
            for k in range(lam_bins):
                if transMatrix[j][k]>0:
                    Q = Q+((rowsumt*((splitMatrices[i][j][k]-transMatrix[j][k])**2))/(transMatrix[j][k]))
    return Q


def findLambda(data,lam_bins=10,numPeriods=10,maxOffset=50):
    """
    Find the time scale (Lambda) in time steps that preserves the Markov
    property.

    Parameters
    ----------
    data : TYPE, numpy array
        Time series being used.
    lam_bins : TYPE, integer.
        Number of bins to identify states used for finding the Markov value.
        The default is 10.
    numPeriods : TYPE, integer
        Number of sections to divide the time series into. The default is 10.
    maxOffset : TYPE, integer.
        The maximum number of time offsets to calculate before stopping.
        The default is 50.

    Returns
    -------
    TYPE integer
        Smallest number of time steps that preserves the Markov condition.

    """

    #find q values over the range of offsets
    Qs = np.zeros(maxOffset)
    for i in range(maxOffset):
        Qs[i] = findQ(data,lam_bins=lam_bins, numPeriods=numPeriods,tau=i+1)
    
    offsets = np.arange(maxOffset)
    Qs = np.log(Qs)
    

    # Find first local max of plotted Q values
    greatest = 0
    count = 0
    smallest = np.argmin(Qs[5:])+5
    greatest=smallest
    for i in range(smallest,len(offsets)):
        if Qs[i]>Qs[greatest]:
            greatest = i
            count = 0
        else:
            count = count+1
        if count==5:
            break
    return (greatest+1)



def KM(TSeries, lambda_1, dt, num_bins = 200, bin_lims = False):
    """
    Calculate the KM coefficients for a time series

    Parameters
    ----------
    TSeries : TYPE, numpy array
        Time series to calculate KM coefficients for.
    lambda_1 : TYPE, integer
        Markov property number of time steps.
    dt : TYPE, float
        Time step size.
    num_bins : TYPE, integer
        Number of bins to bin the data by. The default is 200.
    bin_lims : TYPE, numpy array with len = 2
        Upper and lower bound for bins. The default is False.

    Returns
    -------
    bins : TYPE, 1xnum_bins numpy array
        Bins or x values for KM coefficients.
    D1_e : TYPE, 1xnum_bins numpy array
        D1 coefficients extrapolated to zero.
    D2_e : TYPE, 1xnum_bins numpy array
        D2 coefficients extrapoalted to zero.

    """
    if bin_lims is False:
        bin_lims = [np.min(TSeries)+np.std(TSeries),np.max(TSeries)-np.std(TSeries)]
    if lambda_1 is False:
        lambda_1 = findLambda(TSeries, lam_bins = 10)
    bins = np.linspace(bin_lims[0],bin_lims[1],num_bins)
    dx = np.mean(np.diff(bins))
    bin_center = bins + dx/2
    bins2dx,bins2dy = np.meshgrid(bin_center,bin_center)
    bins_sub = bins2dx - bins2dy
    TSeries_dig = np.digitize(TSeries,bins)
    for tau in range(lambda_1,2*lambda_1):
        m = transition_matrix(TSeries_dig,tau)
        m = m[1:,1:]
        D1_pi = (bins_sub*m)/(tau*dt)
        D2_pi = (bins_sub**2*m)/(2*tau*dt)
        if 'D1' in locals():
            D1 = np.vstack([D1,np.sum(D1_pi,axis=1)/(tau*dt)])
            D2 = np.vstack([D2,np.sum(D2_pi,axis=1)/(tau*dt)])
        if 'D1' not in locals():
            D1 = np.sum(D1_pi, axis=1)/(tau*dt)
            D2 = np.sum(D2_pi, axis=1)/(2*tau*dt)
    
    Tau = np.linspace(lambda_1,2*lambda_1,np.shape(D1)[0])
    D1_e = []
    D2_e = []
    
    # Extract from calculated tau values to zero
    for i in range(np.shape(D1)[1]):
        d1 = D1[:,i]
        if not np.any(d1):
            D1_e = np.append(D1_e,0.0)
        if np.any(d1):
            idx = np.where(d1!=0)
            D1fun = np.poly1d(np.polyfit(Tau[idx],d1[idx],1))
            D1_e = np.append(D1_e,D1fun(0.0))
        d2 = D2[:,i]
        if not np.any(d2):
            D2_e = np.append(D2_e,0.0)
        if np.any(d2):
            idx = np.where(d2!=0)
            D2fun = np.poly1d(np.polyfit(Tau[idx],d2[idx],1))
            D2_e = np.append(D2_e,D2fun(0.0))
    return bins,D1_e,D2_e

def find_KM_fit_coeffs(bins,D1,D2,D1_order=1,D2_order=2):
    """
    Find the polynomial fit used to generate the KM coefficients.

    Parameters
    ----------
    bins : TYPE numpy array
        Center of the bins from KM coefficient calculation (usually considered
        x in written KM equations).
    D1 : TYPE, numpy array
        First-order KM values corresponding to the bins.
    D2 : TYPE, numpy array
        Second-order KM values corresponding to the bins.
    D1_order : TYPE, integer
        Order for fitting function for D1. The default is 1.
    D2_order : TYPE, integer
        Order for fitting function for D2. The default is 2.

    Returns
    -------
    D1_coeffs : numpy array with len = D1_order + 1
        Fitting coefficients for D1 used to generate the time series.
    D2_coeffs : numpy array with len = D2_order + 1
        Fitting coefficients for D2 used to generate the time series.

    """
    idx1 = np.where(D1!=0)
    D1_coeffs = np.polyfit(bins[idx1],D1[idx1],D1_order)
    idx2 = np.where(D2!=0)
    D2_coeffs = np.polyfit(bins[idx2],D2[idx2],D2_order)
    return D1_coeffs, D2_coeffs

def regenerate_ts(x_0,D1_coeffs,D2_coeffs,N=2000,dt=0.1):
    """
    Regenerate a time series from KM coefficients

    Parameters
    ----------
    x_0 : TYPE, float.
        Initial Value.
    D1_coeffs : TYPE, numpy array
        Coefficients for first order polynomial on D1.
    D2_coeffs : TYPE, numpy array
        Coefficients for second order polynomial on D2.
    N : TYPE, TYPE, integer
        How many time steps to solve for. This should be equal to the number
        of time steps in the time series. The default is 2000.
    dt : TYPE, float
        Time step size. The default is 0.1.

    Returns
    -------
    X : Numpy array of dtype: np.float
        Regenerated time series from KM coefficients.

    """
    X_i = x_0
    X = []
    for i in range(N):
        LocalDrift = X_i*D1_coeffs[0]+D1_coeffs[1]
        LocalDiffusion = X_i**2*D2_coeffs[0]+X_i*D2_coeffs[1]+D2_coeffs[2]
        if(LocalDiffusion<0):
            LocalDiffusion = np.abs(LocalDiffusion)
        X_i=X_i+(LocalDrift)*dt+np.sqrt(2.0*LocalDiffusion*dt)*np.random.normal(loc=0,scale=1.0)
        #LocalDiffusions.append(LocalDiffusion)
        #LocalDrifts.append(LocalDrift)
        X = np.append(X,X_i)
    X = np.asarray(X)
    return X
    
