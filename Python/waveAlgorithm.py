#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@brief: Script to calculate the siginificant wave height of a time series of acceleration data from an IMU
@version: 25/09/2019
@author: yusufheylen
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import conditioning as cond
from scipy.stats import pearsonr as corr

from scipy.fftpack import fft, ifft
from scipy.constants import g 
    
def heightTimeSeries(azFiltered, fs, plotFilter=False, plotFFT=False, plotHeight=False ):
    """
        Calculate the height time series from an array of filtered acceleration readings. 
        @param: azFiltered - the acceleration time series to calculate the height of
        @param: fs - the sampling frequency
        @returns: height - the height time series
    """
    df = fs/len(azFiltered)
    fn = fs/2
    
    if (plotFilter or plotFFT ):
        frequencyBinArrayWeights = []
        responseWeightArray = []
    if (plotFFT):
        frequencyBinArray = []

    weightedFFtArray = []
    binNo = 0    
    for fi in fft(azFiltered):
        #Frequency of bin = binNo * ∆f
        weightedFFtArray.append(H(binNo*df, fn)*fi)
        binNo+=1
        
        if (plotFilter):
            if(binNo*df <= 1.0):
                frequencyBinArrayWeights.append(binNo*df)
                responseWeightArray.append(H(binNo*df, fn))
        if (plotFFT):
           frequencyBinArray.append(binNo*df)

    if (plotFilter):
        plt.plot(frequencyBinArrayWeights, responseWeightArray, label="Frequency Weights")
        plt.legend(loc='lower right')
        plt.show()
    if (plotFFT):
        plt.plot(frequencyBinArray, fft(azFiltered), label="FFT")
        plt.legend(loc='lower right')
        plt.show()
     
    timeAxis = []
    height = []
    n = 0
    for e in ifft(weightedFFtArray):
        height.append(e.real)
        timeAxis.append(n * 1/fs)
        n+=1
            
    if(plotHeight):    
        plt.plot(timeAxis,height, label="Wave height")
        plt.legend(loc='lower right')
        plt.show()
        
    return height, timeAxis


def H(f,fn, f1 =  0.02, f2 =  0.03):
    """
        Compute half-cosine taper function.
        @param: fn = Nyquist frequency 
        @param: f = frequency of fourier coefficient
        @param: f1 & f2 = values to restrict effect of transient
    """    
    if ( f < f1 and f > 0):
        return 0
    elif (f <= f2 and f >= f1):
        return 0.5 * ( (1 - np.cos(np.pi * ( (f - f1)/(f2 - f1) ) ) ) * (-1 / ((2*np.pi*f)**2) ) )
    elif ( f > f2 and f < fn):
        return -1/((2*np.pi*f)**2)
    else:
        return 0
    
    
def significantWaveHeight(arrH):
    return 4*np.std(arrH)

def polyfit(t, polycoeff):
    return polycoeff[0] + t*polycoeff[1] + t**2*polycoeff[2]  + t**3*polycoeff[3] + t**4*polycoeff[4] + t**5*polycoeff[5] + t**6*polycoeff[6] + t**7*polycoeff[7] + t**8*polycoeff[8] + t**9*polycoeff[9] + t**10*polycoeff[10] + t**11*polycoeff[11] + t**12*polycoeff[12] + t**13*polycoeff[13] +t**14*polycoeff[14] + t**15*polycoeff[15] + t**16*polycoeff[16] + t**17*polycoeff[17] + t**18*polycoeff[18] + t**19*polycoeff[19]

def spectralSignificantWaveHeight(heights, fs):
    f, pxx = sig.periodogram(heights, fs)
    s = 0
    for i in range(1, 1025):
        s+=pxx[i]
    s*= fs/len(heights)
    return 4*np.sqrt(s)
    
    
def main():


    
    #To plot yost v theoretical
#    arrPitch = []
#    arrRoll = []
#    arrAx = []
#    arrAy = []
    
#    imu = int(input("Enter 0 for YOST IMU, 1 for MPU6050:\n"))
#    fileName = input("Enter the acceleration reading file path and name:\n")
    
    
    experiments = ["../Data/YOST_stewart_0degPitch_10sPeriod_test_1.txt","../Data/MPU6050_stewart_0degPitch_10sPeriod_test_1.txt",
                   "../Data/YOST_stewart_0degPitch_20sPeriod_test_2.txt","../Data/MPU6050_stewart_0degPitch_20Period_test_2.txt",
                   "../Data/YOST_stewart_20degPitch_20sPeriod_test_3.txt","../Data/MPU6050_stewart_20degPitch_20Period_test_3.txt"]
    plot = True
    displacements = []
    sigWaveHeights = []
    
    for i in range(0,6) :
        arrAz = []
        totalTime = 0
       
        imu = i%2   #YOST => 0, MPU6050 => 1
        with open(experiments[i]) as f:
            #Valid files
            #YOST_stewart_0degPitch_10sPeriod_test_1.txt 
            #YOST_stewart_0degPitch_20sPeriod_test_2.txt
            #YOST_stewart_20degPitch_20sPeriod_test_3.txt
            #MPU6050_stewart_0degPitch_10sPeriod_test_1.txt
            #MPU6050_stewart_0degPitch_20Period_test_2.txt
            #MPU6050_stewart_20degPitch_20Period_test_3.txt
            
            #If YOST IMU (imu = 0)
            #Data format: "%int(Month)/%int(Day)/%int(Year),%int(Hours):%int(Minutes):%float(Seconds),
            # %float(OrientPitch),%float(OrientYaw),%float(OrientRoll),
            # %float(CorrectedGyroX),%float(CorrectedGyroY),%float(CorrectedGyroZ),
            # %float(CorrectedAccelX),%float(CorrectedAccelY),%float(CorrectedAccelZ),
            # %float(CorrectedMagX),%float(CorrectedMagY),%float(CorrectedMagZ)"
            if (imu == 0):
                f.readline() # Read in first line - this is the Foramt 
                
                #Get values from file
                startTime = 0 
                endTime = 0
                for line in f:  
                    
                    row = line.split(',')
                    
                    #Get start time
                    if(startTime == 0):
                        startTime = row[0].split(' ')[1] 
                    
                    #Get end time
                    endTime = row[0].split(' ')[1]
                    
        
                    #Select relevent accleration data - comment out if plotting yost v theoretical
                    row = row[7:10]
    
        
                    #Set upper bound of 0.5g Az
                    if (float(row[1]) > 0.5*g):
                        row[1] = str(0.5*-g)
                    arrAz.append(float(row[1])*-g ) #comment out if plotting yost v theoretical

                    #This is also used to compare yost with the true signal
#                    arrAz.append(float(row[-5])*-g )
#                    arrAx.append(float(row[-6])*-g )
#                    arrAy.append(float(row[-4])*-g )   
#                    arrPitch.append(float(row[1]))
#                    arrRoll.append( float(row[3]) )
                
                #Calculate the sampling frequency
                startTime = startTime.split(':')
                endTime = endTime.split(':')
                totalTime = []
                totalTime.append(float(endTime[0]) - float(startTime[0]))
                totalTime.append(float(endTime[1]) - float(startTime[1]))
                totalTime.append(float(endTime[2]) - float(startTime[2]))
                
                totalTime = totalTime[0]*60*60 + totalTime[1]*60 + totalTime[2]
                
            #Else MPU6050 (imu = 1)
            #Data format: "int(timeSinceStart ms), float(accelAx mg), float(accelAy mg), float(accelAz g)"
            else:
                startTime = -1
                endTime = 0
                for line in f:
                    #Format is: int ms, float ax, float ay, float az
                    row = line.split(',')
                    if(startTime == -1 ):
                        startTime = float(row[0])*10**-3
                    endTime = float(row[0])*10**-3
                    #Set upper bound of 0.5g Az
                    if (float(row[3]) > 0.5*g):
                        row[1] = str(0.5*-g)
                    #arrAx.append(float(row[1])*-g/1000 )
                    #arrAy.append(float(row[2])*-g/1000 )   
                    arrAz.append(float(row[3])*-g )
    
                totalTime = endTime - startTime
                    
                        
            
            fs = len(arrAz)/(totalTime) #Sampling frequency
            fs = round(fs)  #Account for errors
            
            ##Debuging and graphing
            #print("Sampling rate = " + str(fs))
            #trueVerticalAcceleration(arrAx, arrAy, arrAz, arrPitch,arrRoll, fs)
            ##EndDebug
            
        
            #Condition signal:
            azFiltered = cond.condition(arrAz, fs, plot)
            
            #Calculate Wave height time series
            eta, times = heightTimeSeries(azFiltered, fs, plot, plot, plot)        
           
            #Resample to allow for comparison between the imus (has to have same amount of samples)
            eta180, times = sig.resample(eta,180,  t=times)
            if (plot):
                plt.plot(times, eta180, label="Reasmpled heights")
                plt.legend(loc='lower right')
                plt.show()
    
            displacements.append(eta180)
            
            ht = significantWaveHeight(eta)
            hs = spectralSignificantWaveHeight(eta, fs)
            sigWaveHeights.append( (ht,hs) )
            
            
#    print(displacements)
    h = 0.045
    c = 0.155
    f = 0.1 
    t = np.arange(0,90,0.5)
    s = h*np.sin(2*np.pi*f*t);

    for j in range(0,6):
        if (j%2 == 0):
            print("YOST Significant Wave Height (Ht, Hs) for test " + str(round(j*2/5)) +  ": Ht=" + '{:6f}'.format(sigWaveHeights[j][0]*1000) + "mm Hs=" + '{:6f}'.format(sigWaveHeights[j][1]*1000) )
        else:
            print("MPU6050 Significant Wave Height(Ht, Hs) for test " + str(round(j*2/5)) + ": Ht=" + '{:6f}'.format(sigWaveHeights[j][0]*1000) +  "mm Hs=" + '{:6f}'.format(sigWaveHeights[j][1]*1000) )
        
    print("Theoretical Significant Wave Height: " + '{:6f}'.format(significantWaveHeight(s)*1000) + "mm" )
    
    for k in range(0,6,2):
        print("Pearson coerrelation coefficient between IMUs for test " + str(int(k/2)) + " is: " + '{:6f}'.format(abs(corr(displacements[k],displacements[k+1])[0])))
                  

if __name__ == "__main__":
    main()
    
   
#######################################################################################################################
#
##    This is used to compare yost imu to theoretical wave    
#    
#def trueVerticalAcceleration(arrAx, arrAy, arrAz, arrPitch,arrRoll, fs, plot=True ):
#    tva = []
#    for i in range(len(arrAx)):
#        tva.append(cond.calculateTrueVerticalAcceleration(arrAx[i], arrAy[i], arrAz[i],arrPitch[i], arrRoll[i]))
#    tvaFiltered  = cond.condition(tva, fs, decimate=True)
#    eta, times = heightTimeSeries(tvaFiltered, fs)
#    eta180, a = sig.resample(eta,180,  t=times)
#    if(plot):
#        t = np.arange(0,90,0.5)
#        h = 0.045
#        f = 0.1 
#        s = h*np.sin(2*np.pi*f*t);
#        plt.plot(t, s, 'b', label="Theoretical wave displacement")
#        plt.plot(t, eta180, 'r', label="Calcultated Wave displacement")
#        plt.legend(loc='lower right')
#    plt.show()
#    print("Ht " + str(significantWaveHeight(eta)) )
#    print("Hs " + str(spectralSignificantWaveHeight(eta, fs)) )
    