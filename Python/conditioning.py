#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@brief: Script to condition a time series of acceleration data from an IMU for ocean wave analysis
@version: 25/09/2019
@author: yusufheylen
"""
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

from copy import deepcopy
from fractions import Fraction

def removeSpikes(arr):
    """
        Remove spikes in data that is above/below six stds from the mean
        @param: arr array to remove spikes from
        @returns: arr = array with removed spikes and (linear) interpolated values
    """
    for i in range (0,3):
        spikes = 0
        mean = np.mean(arr)
        stdDiv = np.std(arr)
        i = 0
        spikeIndexArr = []
        for reading in arr:
            if (reading > (mean + stdDiv*3) or reading < (mean - 3*stdDiv)):
                spikes += 1
                spikeIndexArr.append(i)
            i += 1
        #Interpolate
        for i in spikeIndexArr:
            arr[i] = (arr[i-1] + arr[i+1])/2
    return arr

def detrend(arr, fs):
    """
        Remove trends in time series. Equivalant to a high pass filter
        with time constant = (1/fs)/(1 - k)
        @param: arr = array to remove trends from
        @param: fs = sampling frequency
        @returns: detArr = detrended array
    """
    mean = np.mean(arr)
    detArr = []
    sn_prev = 0
    k = 1 - (1/fs)/100
    for yn in arr:
        sn = (yn-mean) + k*sn_prev
        detArr.append((yn-mean) - (1-k)*sn)
        sn_prev = sn
    return detArr
    

def analogLPF(arr, fn, fc = 8, n=1):
    """
        Apply a 1st order digital low pass filter on an array
        @param: arr = signal to filter
        @param: fn = nyquist frequency
        @param: fc = cut-off frequency
        @returns: filtered signal        
    """
    w = fc/fn
    b, a = sig.butter(n,w, 'low', True)
    return sig.lfilter(b,a, arr)
    
    
def decimate(arrAz, fs):
    """
        Decimate the signal to 2Hz
        @param: arrAz - the signal to decimate
        @param: fs - the sampling frequency
        @param: target - the frequency wanted
        @returns: decimated array
    """

    #upsample to 640Hz
    frac = Fraction(640/fs).limit_denominator()
    s640 = sig.resample_poly(arrAz, frac.numerator, frac.denominator )
    s640 = 9*analogLPF(s640, fs/2)

    #Decimate to 80Hz
    frac = Fraction(80/fs).limit_denominator()
    s80 = sig.resample_poly(s640, frac.numerator, frac.denominator )
    s80 = 9*analogLPF(s80, 80/2, 1, 2)
#    
#    #Decimate to 2Hz
    frac = Fraction(2/fs).limit_denominator()
    s2 = sig.resample_poly(s80, frac.numerator, frac.denominator )
    s2 = 9*analogLPF(s2, 2/2, 1, 2)
    
    return s2

    

def digitalHPF(arr, fn, fc = 0.025):
    """
        Digital high pass filter in time domain to be used as alternative to half-cosine taper in frequency domain,
        for debugging and comparison purposes
        @param: arr = signal to filter
        @param: fn = nyquist frequency 
        @param: fc = cut-off frequency
        @returns: Filtered signal (array)
    """
    w= fc/fn
    b, a = sig.butter(1, w, 'hp')
    return sig.lfilter(b,a,arr)
    
    

def condition(accelArr, fs, plot=False, dec=False):
    """
        Apply conditioning to signal and return conditioned signal as well as the times.
        @param: accelArr - acceleration time series to condition
        @param: fs - sampling frequency
        @param; plot - specify if to plot
        @param: dec - specify if to upsample and decimate
        @returns: accelArr conditioned, timeAccelAxis corresponding times
    """    
    if (plot == True):
        rawAz = deepcopy(accelArr)
    
    accelArr = removeSpikes(accelArr)
    accelArr = detrend(accelArr,fs)


    

   
    timeAccelAxis = []
    t = 0
    for i in range(len(accelArr)):
        timeAccelAxis.append(t * 1/fs)
        t +=1
        
    if(plot == True):
        plt.plot(timeAccelAxis, rawAz, label="Raw accel.")
        plt.legend(loc='lower right')
        plt.show()
        plt.plot(timeAccelAxis, accelArr, label="Conditioned accel.")
        plt.legend(loc='lower right')
        plt.show()

    
    if(not(dec)):
        accelArr = analogLPF(accelArr, fs/2, fc=1)
    
        if(plot == True):
            plt.plot(timeAccelAxis, accelArr, label="Filtered accel.")
            plt.legend(loc='lower right')
            plt.show()
            
    if(dec):
        accelArr = decimate(accelArr, fs)
    
    return accelArr


#######################################################################################################################
#
##    THIS IS A POSSIBLE EXTENSION IF A GYRO SEM HAS BEEN DEVELOPED
##    NOTE: THE calculateTrueVerticalAcceleration FUNCTION WOULD BE CALLED *BEFORE* APPLYING LPF
#        
#def calculateTrueVerticalAcceleration(Ax, Ay, Az, pitch, roll):
#    """
#        Compute the true vertical acceleration relative to the earth. If tilt < 10˚  => use Az only. 
#        Else calculate tva. TODO REPLACE ACCEL TO ANGLE VS GYRO READING 
#        @param: A* = Acceleration mesured in the * direction
#        @param: theta = pitch reading at same timestamp
#        @param: phi = roll reading at same timestap
#        
#    """
#    if(pitch > 10*(np.pi/180) or roll > 10*(np.pi/180)):
#        return (np.sin(pitch)*Ax 
#                + np.sin(roll)*np.cos(pitch)*Ay 
#                + np.cos(roll)*np.cos(pitch)*Az
#        )
#    else:
#        return Az
#
#
#   NOTE: THE TILT FUNCTION WOULD HAVE TO BE REDEFINED TO CALCULATE TILT FROM THE PITCH AND ROLL         
#def calculateTilt(Gx, Gy, Gz):
#    """
#        Calculate the tilt - rho
#    """
#    return np.arccos(Gz / (np.sqrt(Gx**2 + Gy**2 + Gz**2)))
#
##   NOTE: THIS FUNCTION WOULD BE DEPRECIATED / REMOVED 
#def calculatePitch(Gx, Gy, Gz):
#    """
#        Calculate the pitch - theta
#    """
#    return np.arctan( -Gx / ( np.sqrt( Gy**2 + Gz**2 ) ) )
# 
##   NOTE: THIS FUNCTION WOULD BE DEPRECIATED / REMOVED 
#def calculateRoll(Gy, Gz):
#    """
#        Calculate the roll - phi
#    """
#    return np.arctan(Gy / Gz)
#        
#######################################################################################################################