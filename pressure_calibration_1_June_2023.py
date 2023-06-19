#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 16:33:13 2023

@author: Mohammad-Amin Aminian

DFG Calibration

"""
from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import coincidence_trigger
from obspy.signal.trigger import classic_sta_lta
import matplotlib.pyplot as plt
import scipy
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import numpy as np
import obspy
import tiskit
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import coincidence_trigger
from obspy.signal.trigger import plot_trigger

def calculate_spectral_ratio(stream,mag = 7,freq1=0.05,freq2=0.001,rho=1028,nseg = 2**10,TP = 5):
    '''
    Parameters
    ----------
    stream : TYPE
        Stream Before Removeing instrument response.
    mag : TYPE, optional
        Minimumm Magnitude for earthquakes to use in calulation. The default is 7.
    freq1 : TYPE, optional
        high frequency to filter the data . The default is 0.05.
    freq2 : TYPE, optional
        low frequency to filter the data. The default is 0.001.
    rho : TYPE, optional
        density of ocean water. The default is 1028.
    nseg : TYPE, optional
        number of segments to calculate coherence funtion. The default is 2**10.
    TP : TYPE, optional
        Length of taper at one end in minutes. The default is 5.

    Returns
    -------
    ratio : TYPE
        DESCRIPTION.

    '''
    stream2 = stream.copy()
    stream2.clear()
    
    invz = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BHZ",
        location="*",
        level="response")

    invp = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BDH",
        location="*",
        level="response")
    eq_spans = tiskit.TimeSpans.from_eqs(stream.select(channel='*Z')[0].stats.starttime,
                                     stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=mag, days_per_magnitude=0.5)
    
    for i in range(0,len(eq_spans)):
        stream1= stream.copy()
        stream2 =stream2 + stream1.trim(eq_spans.start_times[i-1],eq_spans.start_times[i-1]+2*3600)
    
    stream2.sort(['starttime','channel'])
    # stream.trim(eq_spans.start_times[i-1]+2200,eq_spans.start_times[i-1]+4000)
    
    stream2.select(channel="*Z").remove_response(inventory=invz,
                                            output="VEL", plot=False)

    stream2.select(channel="*H").remove_response(inventory=invp,
                                            output="DEF", plot=False)
    stream2.filter("lowpass", freq=freq1)
    stream2.filter("highpass", freq=freq2)

    stream2.taper(max_percentage=0.1,type='hann',max_length=(TP*60))
    stream2.sort(['starttime','channel'])

    spectrum_z = np.zeros([len(eq_spans),stream2[0].stats.npts])
    spectrum_p = np.zeros([len(eq_spans),stream2[0].stats.npts])
    ratio = np.zeros([len(eq_spans),stream2[0].stats.npts])
    Czp = np.zeros([len(eq_spans),int(nseg/2 + 1)])
    acc_p = np.zeros([len(eq_spans),stream2[0].stats.npts])
    freqs = np.fft.fftfreq(len(stream2.select(channel="*Z")[0].data),d=stream2.select(channel="*Z")[0].stats.sampling_rate)

    for i in range(0,len(eq_spans)):
        f,Czp[i] = scipy.signal.coherence(stream2.select(component='Z')[i].data,
                                         stream2.select(component='H')[i].data,
                                         fs=stream2[i].stats.sampling_rate,
                                         nperseg =nseg,noverlap=(nseg*0.5),
                                         window=scipy.signal.windows.tukey(nseg,
                                         (TP*60*stream2[i].stats.sampling_rate)/nseg))
    
        acc_p[i] =  stream2.select(channel="*H")[i].data / (rho * -invz[0][0][0].elevation) 
        # in code this is velocity , I'll convert it when im computing spectrum_p

    # Perform FFT on the signals
        spectrum_z[i] = np.fft.fft(stream2.select(channel="*Z")[i].data)* (np.pi *2*freqs)
        spectrum_p[i] = np.fft.fft(stream2.select(channel="*H")[i].data)

        spectrum_p[i] = np.fft.fft(acc_p[i]) 
        ratio[i] = np.abs(spectrum_p[i]) / np.abs(spectrum_z[i])
        
    # Calculate the spectral ratio
    
    for i in range(0,len(eq_spans)):
        plt.rcParams.update({'font.size': 20})
        plt.figure(dpi=300,figsize=(10,15))
        
        # Plot 2: Acceleration
        plt.subplot(411)
        plt.plot(stream2.select(channel="*Z")[i].times(), stream2.select(channel="*Z")[i].data,linewidth=3,color='blue')
        plt.xlabel('Time (s)')
        plt.ylabel('Vertical Acc[m/s^2]')
        plt.grid(True)
        # Plot 2: Pressure
        plt.subplot(412)
        plt.plot(stream2.select(channel="*H")[i].times(), stream2.select(channel="*H")[i].data,linewidth=3,color='blue')
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure [pa] ')
        # plt.title('Stream2')
        plt.grid(True)


    # Plot 3: Coherence
        plt.subplot(413)
        plt.semilogx(f,Czp[i],linewidth=3,color='blue')
        plt.ylabel("Coherence ")
        plt.xlabel("Frequency ")
        plt.grid(True)
        plt.ylim([0,1])
        plt.xlim([0.001,freq1*2])
        plt.grid(True)

    # Calculate spectral ratio

    # Plot 4: Spectral Ratio/Frequency
        plt.subplot(414)
        positive_freq_indices = np.where((freqs) >= 0)

    # Plot the right-sided spectrum
        plt.loglog(freqs[positive_freq_indices],(ratio[i][positive_freq_indices]),linewidth=2,color='blue')
        plt.xlim([0.001,freq1*2])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Spectral Ratio')
        plt.title('Spectral Ratio ')
        plt.grid(True)
    
        plt.tight_layout()
        plt.show()
    
    # Calculate the spectral ratio
    
    plt.rcParams.update({'font.size': 20})
    plt.figure(dpi=300,figsize=(12,15))
    plt.subplot(411)
    plt.plot(stream.select(channel="*Z")[0].times(), stream.select(channel="*Z")[0].data,linewidth=3)
    plt.xlabel('Time (s)')
    plt.ylabel('Vertical Acc[m/s^2]')

    
    # Plot 2: Stream2
    plt.subplot(412)
    plt.plot(stream.select(channel="*H")[0].times(), stream.select(channel="*H")[0].data,linewidth=3)
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure [pa] ')
    # plt.title('Stream2')


    # Plot 3: Stream2 as acceleration
    plt.subplot(413)
    for i in range(0,len(eq_spans)):
        plt.semilogx(f,Czp[i],linewidth=3)
    plt.ylabel("Coherence ")
    plt.xlabel("Frequency ")
    plt.grid(True)
    plt.ylim([0,1])
    plt.xlim([0,freq1*2])

    # Calculate spectral ratio

    # Plot 4: Spectral Ratio/Frequency
    plt.subplot(414)
    positive_freq_indices = np.where((freqs) >= 0)

    # Plot the right-sided spectrum
    for i in range(0,len(eq_spans)):
        if np.median(Czp[i]) > 0.4 :
            plt.semilogx(freqs[positive_freq_indices],(ratio[i][positive_freq_indices]),linewidth=0.2,color='r')
            print(i)
            # ratio2 = ratio2 + ratio[i]
            
    plt.semilogx(freqs[positive_freq_indices],np.mean(ratio[:],axis=0)[positive_freq_indices],linewidth=1,color='blue')
    # plt.semilogx(freqs[positive_freq_indices],ratio2[positive_freq_indices]/3,linewidth=1,color='blue')
    plt.xlim([freq2,freq2*3])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Spectral Ratio')
    plt.title('Spectral Ratio ')
    # plt.ylim([0,100])
    
    plt.tight_layout()
    plt.show()
    
    plt.figure(dpi=300)
    plt.plot(stream2[0].times(),stream2[3].data/np.max(stream2[3].data))
    plt.plot(stream2[0].times(),stream2[0].data/np.max(stream2[0].data))
    # plt.xlim([9.7*10e4,9.8*10e4])
    plt.xlabel('Time[s]')
    plt.ylabel('Normilized Acc or Pressure')

    
    plt.figure(dpi=300)
    plt.specgram(stream[3].data,Fs = 2.1,NFFT=256)
    plt.ylim([0.001,0.05])
    plt.xlabel('Time[s]')
    plt.ylabel('Frequency[Hz]')

    return ratio

#%%
# Calibrating via P-wave Arrival
client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR38"



def pressure_calibration(stream,mag=7,i=1):
    
    invz = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BHZ",
        location="*",
        level="response")

    invp = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BDH",
        location="*",
        level="response")
    eq_spans = tiskit.TimeSpans.from_eqs(stream.select(channel='*Z')[0].stats.starttime,
                                     stream.select(
                                         channel='*Z')[0].stats.endtime,
                                     minmag=mag, days_per_magnitude=0.5)
    
    stream.trim(eq_spans.start_times[i-1],eq_spans.start_times[i-1]+10*3600)
    
    stream.select(channel="*Z").remove_response(inventory=invz,
                                            output="ACC", plot=False)

    stream.select(channel="*H").remove_response(inventory=invp,
                                            output="DEF", plot=False)


    df = stream[0].stats.sampling_rate

    cft = obspy.signal.trigger.recursive_sta_lta(
        stream.select(channel="*Z")[0].data, int(5 * df), int(10 * df))
    
    # trig = coincidence_trigger("recstalta", 7, 2, stream, 1, sta=int(5 * df), lta=int(10 * df))
    
    p_arrival_times = obspy.signal.trigger.trigger_onset(cft,1,0.2) / df   
    plot_trigger(stream.select(channel="*Z")[0], cft,1, 0.5)

    
    stream.trim(stream[0].stats.starttime + p_arrival_times[i-1][0] -2 ,stream[0].stats.starttime + p_arrival_times[i-1][0]+2)
    stream.detrend()
    rho = 1028
    Acc = stream[0].data / (rho*-invz[0][0][0].elevation)
    
    plt.plot(stream[3].data),plt.plot(Acc)
    
    plt.plot(stream[3].data / -Acc),plt.plot(stream[3].data)
    
    scipy.stats.pearsonr(stream[0].data,stream[3].data)
