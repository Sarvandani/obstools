#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:59:48 2023

@author: mohammadamin
Downsampling OBS data

"""

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import tiskit
import scipy
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import obspy
import obspy.signal
nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()



client = Client("RESIF")
net = "YV"
sta = "RR38"
t1 = UTCDateTime("2012-12-15T12:00:00.000000Z")
t2 = t1 + 24 * 3600


inv = client.get_stations(
    network=net,
    station=sta,
    channel="*",
    location="*",
    level="response")


stream = client.get_waveforms(net, sta, "*", "*", starttime= t1,endtime=t2)

decim = tiskit.Decimator([5, 3, 2])
stream_decim = decim.decimate(stream)
inv_decim = decim.update_inventory(inv, stream)
print(inv)

stream.remove_response(inventory=inv,output="DISP", plot=False)

stream_decim.remove_response(inventory=inv_decim,output="DISP", plot=False)

compare_z = stream.select(channel='*Z') + stream_decim.select(channel='*Z')

compare_z.plot(dpi = 300 ,method='full',equal_scale=True,figsize=(12,8))


nseg = 2**17
TP = 5

f1,Dz = scipy.signal.welch(stream.select(component='Z')[0],
                              fs=stream[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg*0.9),
                              window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg))

nseg = 2**12
f,Dz_d = scipy.signal.welch(stream_decim.select(component='Z')[0],
                              fs=stream_decim[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg*0.9),
                              window=scipy.signal.windows.tukey(nseg,(TP*60*stream_decim[0].stats.sampling_rate)/nseg),nfft=2**17)

plt.rcParams.update({'font.size': 20})
plt.figure(dpi=300,figsize=(12,8))
plt.plot(f1,10*np.log10(Dz* (2*np.pi*f1**2)),label="BHZ",linewidth=3,color='r')             
plt.plot(f,10*np.log10(Dz_d* (2*np.pi*f**2)),label="BHZ Decimated",linewidth=2,color='g')             

plt.xscale('log')
plt.title(stream[0].stats.network+"."+stream[0].stats.station+"  "+str(stream[0].stats.starttime)[0:10]+"--"+str(stream[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acceleration$((m/s^2)^2/Hz)$ [dB] ")
plt.grid(True)
plt.ylim([-200 ,-80])
plt.xlim([0.001,1])
plt.tight_layout()

plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=12)

#%%
from tiskit.rptransient import  Transients, PeriodicTransient as PT

Rot = tiskit.CleanRotator(stream_decim, remove_eq=True,filt_band=(0.001, 0.1))
    
stream_decim_rot = Rot.apply(stream_decim,horiz_too=True)
    
    
transients = {'RR28': [PT("1h", 3620.3, 0.05, [-820, 330], stream[0].stats.starttime)],
                  'RR29': [PT("1h", 3620.3, 0.05, [-290, 135], stream[0].stats.starttime)],
                  'RR38': [PT("1h", 3620.3, 0.05, [-820, 330], stream[0].stats.starttime)]}

rt = Transients(transients[sta])
    
zdata = stream_decim_rot.select(channel='*Z')[0]
       
zdata.detrend(type='demean')
zdata.detrend(type='linear')
zdata.detrend(type='demean')
 
eq_spans = tiskit.TimeSpans.from_eqs(zdata.stats.starttime, zdata.stats.endtime, minmag=5.5,days_per_magnitude=1.5)
    
rt.calc_timing(zdata, eq_spans)
    
rt.calc_transients(zdata, eq_spans, plot=False)
    
cleaned = rt.remove_transients(zdata, plot=False, match=False, prep_filter=True)
    
stream[3] = cleaned
    
zdata.detrend(type='demean')
zdata.detrend(type='linear')
zdata.detrend(type='demean')
    
stream = eq_spans.zero(stream)
    
stream.select(channel='*H')[0].stats.channel = 'BDH'
stream.select(channel='*1')[0].stats.channel = 'BH1'
stream.select(channel='*2')[0].stats.channel = 'BH2'
stream.select(channel='*Z')[0].stats.channel = 'BHZ'
    
#%%
client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR29"

invz = client.get_stations(
    network=net,
    station=sta,
    channel="BHZ",
    location="*",
    level="response")

import os
from obspy import read

B = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/')
B.sort()

A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/')
A.sort()

i = 3

stream = read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+A[i])

i = 2
stream_Raw = read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Decimated/'+B[i])
stream
stream_Raw

stream.trim(t1,t2)
stream_Raw.trim(t1,t2)

# stream.plot(equal_scale=False)

stream.select(channel="*Z").remove_response(inventory=invz,
                                            output="DISP", plot=False)


stream_Raw.select(channel="*Z").remove_response(inventory=invz,
                                            output="DISP", plot=False)


nseg = 2**12
f1,Dz = scipy.signal.welch(stream.select(component='Z')[0],
                              fs=stream[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg*0.9),
                              window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg))

nseg = 2**12
f,Dz_d = scipy.signal.welch(stream_Raw.select(component='Z')[0],
                              fs=stream[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg*0.9),
                              window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg))

plt.rcParams.update({'font.size': 20})
plt.figure(dpi=300,figsize=(12,8))
plt.plot(f1,10*np.log10(Dz* (2*np.pi*f**2)),label="Glitches Removed",linewidth=3,color='g')             
plt.plot(f,10*np.log10(Dz_d * (2*np.pi*f**2)),label="Raw",linewidth=2,color='r')             

plt.xscale('log')
plt.title(stream[0].stats.network+"."+stream[0].stats.station+"  "+str(stream[0].stats.starttime)[0:10]+"--"+str(stream[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acceleration$((m/s^2)^2/Hz)$ [dB] ")
plt.grid(True)
plt.ylim([-200 ,-80])
plt.xlim([0.001,1])
plt.tight_layout()

plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=12)

    

compare_z = stream.select(channel='*Z') + stream_Raw.select(channel='*Z')

compare_z.detrend(type='demean')
compare_z.detrend(type='linear')
compare_z.detrend(type='demean')

compare_z[1].stats.channel = "BHZ-T"

compare_z.plot(dpi = 300 ,method='full',equal_scale=True,figsize=(12,8))


#%%
client = Client("RESIF")
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR38"

invz = client.get_stations(
    network=net,
    station=sta,
    channel="BHZ",
    location="*",
    level="response")

import os
from obspy import read

A = os.listdir("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/')
A.sort()

i = 3

stream = read("/Users/mohammadamin/Desktop/Data/YV/"+sta+'/Transients/'+A[i])



# stream.trim(t1,t2)

# stream.plot(equal_scale=False)

stream.select(channel="*Z").remove_response(inventory=invz,
                                            output="ACC", plot=False)

eq_spans = tiskit.TimeSpans.from_eqs(stream[0].stats.starttime, stream[0].stats.endtime, minmag=5.8,days_per_magnitude=1.5)

stream_EQ = eq_spans.zero(stream)
    
compare_z = stream.select(channel='*Z') + stream_EQ.select(channel='*Z')

compare_z.plot(dpi = 300 ,method='full',equal_scale=True,figsize=(10,8))

nseg = 2**12
f1,Dz = scipy.signal.welch(stream.select(component='Z')[0],
                              fs=stream[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg*0.9),
                              window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg))

nseg = 2**12

data = compare_z[1].data
non_zero_data = data[data != 0]  # Exclude zero values from the data

f,Dz_d = scipy.signal.welch(non_zero_data,
                              fs=stream[0].stats.sampling_rate,
                              nperseg = nseg, 
                              noverlap=(nseg*0.9),
                              window=scipy.signal.windows.tukey(nseg,(TP*60*stream[0].stats.sampling_rate)/nseg))


plt.rcParams.update({'font.size': 20})
plt.figure(dpi=300,figsize=(12,8))
plt.plot(f1,10*np.log10(Dz* (2*np.pi*f**2)),label="BHZ",linewidth=3,color='r')             
plt.plot(f,10*np.log10(Dz_d * (2*np.pi*f**2)),label="BHZ - Earthquakes",linewidth=2,color='g')             

plt.xscale('log')
plt.title(stream[0].stats.network+"."+stream[0].stats.station+"  "+str(stream[0].stats.starttime)[0:10]+"--"+str(stream[0].stats.endtime)[0:10])    
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acceleration$((m/s^2)^2/Hz)$ [dB] ")
plt.grid(True)
plt.ylim([-200 ,-80])
plt.xlim([0.001,1])
plt.tight_layout()
plt.plot(1/nhnm[0],nhnm[1],'--k',label = "New Low/High Noise Model" )
plt.plot(1/nlnm[0],nlnm[1],'--k')
plt.legend(loc = "upper right",fontsize=12)

#%%STA/LTA
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import coincidence_trigger
from obspy.signal.trigger import plot_trigger

df = stream[3].stats.sampling_rate

# cft = classic_sta_lta(stream[3].data, int(2 * df), int(50 * df))

cft = obspy.signal.trigger.recursive_sta_lta(
    stream[3].data, int(10 * df), int(200 * df))


plt.figure(dpi=300)
plot_trigger(stream[3], cft, 5, 4)

# stream.trigger("recstalta", sta=20, lta=10)

# trig = coincidence_trigger("recstalta", 5, 1, stream,
#                            0.1, sta=4, lta=100, details=True)

# trig[0]['time']




