import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import stats
import scipy as sp
from operator import itemgetter
from itertools import groupby
#from sklearn.feature_selection import VarianceThreshold
#from sklearn import cluster
#from sklearn import metrics
#from sklearn.datasets.samples_generator import make_blobs
#from sklearn.preprocessing import StandardScaler
import heapq
#import cPickle
import os


def raster(class_index,ts,window,step):
    ## step in data points 
    ## window: 1s 20000 data points
    ## step: 2ms 40 data points
    raster_index=[]
    blocks=window/step
    for time in ts:
        raster_index_dummy=[]
        for ii in range(blocks):
            step_block=range((int(time)+ii*step),(int(time)+(ii+1)*step))
            for index in class_index:
                if index in step_block:
                    raster_index_dummy.append(ii)
        
        raster_index.append(raster_index_dummy)   
    raster_index=np.asarray(raster_index)       
    return raster_index ##2d arrary of event index
                                        
def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    y = signal.filtfilt(b, a, data)
    return y
   
def extrac_peak_index(data,thr_low,thr_high):
    slope = data[1:] - data[:-1]
    indices = [i+1 for i in range(len(slope)-1) if (slope[i] > 0 and slope[i+1] < 0)]
    extracted_peak_index=[]
    for ele in indices:
        if ((data[ele]>thr_low) and (data[ele]<thr_high)):
            extracted_peak_index.append(ele)
    return extracted_peak_index ## list
 
#def extrac_peak_index(data,thr_low,thr_high):
#    slope = data[1:] - data[:-1]
#    indices = [i+1 for i in range(len(slope)-1) if ((slope[i] > 0 and slope[i+1] < 0) or (slope[i] < 0 and slope[i+1] > 0))]
#    extracted_peak_index=[]
#    for ele in indices:
#        if (((data[ele]>thr_low) and (data[ele]<thr_high)) or ((data[ele]<-thr_low) and (data[ele]>-thr_high))):
#            extracted_peak_index.append(ele)
#    return extracted_peak_index ## list

#define DWT using Haar wavelets 
def discreteHaarWaveletTransform(x):
    N = len(x)
    output = [0.0]*N

    length = N >> 1
    while True:
        for i in xrange(0,length):
            summ = x[i * 2] + x[i * 2 + 1]
            difference = x[i * 2] - x[i * 2 + 1]
            output[i] = summ
            output[length + i] = difference

        if length == 1:
            return output

        #Swap arrays to do next iteration
        #System.arraycopy(output, 0, x, 0, length << 1)
        x = output[:length << 1]

        length >>= 1

## trim the spike trace matrix, remove any that did not peak at 10
def spike_trace_selection(data,index):
    new_data=[]
    new_index=[]
    for ii in range(len(index)):
        if data[ii].max()==data[ii][10]:
            new_data.append(data[ii])
            new_index.append(index[ii])
    new_data=np.asarray(new_data)
    new_index=np.asarray(new_index)
    return new_data,new_index
    
## select trial index with max and go and select trial with nogo
def trial_selection(trial_info):
    go_index=[]
    nogo_index=[]
    for ii in range(trial_info.shape[0]):
        if ((trial_info[ii,0].astype(float)==trial_info[:,0].astype(float).max()) and (trial_info[ii,2]=='GO_REMIND' or trial_info[ii,2]=='GO')):
            go_index.append(ii)
        if (trial_info[ii,2]=='NOGO' or trial_info[ii,2]=='NOGO_REPEAT'):
            nogo_index.append(ii)
    go_index=np.asarray(go_index)
    nogo_index==np.asarray(nogo_index)
    return go_index,nogo_index
    
## function for calcualte accuracy
def predict_accuracy (answer, output):
    
    length = len(answer)
    
    count = 0.0;
    for ii in range(0,length):
        if answer[ii] == output[ii]:
            count = count + 1
    
    return count/length
    
## function for making wavelet transform matrix
            
        
    
        



        
        
    