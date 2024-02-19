import matplotlib.pyplot as plt
import numpy as num

from pyrocko import util, model, io, trace, moment_tensor, gmtpy
from pyrocko import pz
from pyrocko import orthodrome as od
from pyrocko.io import quakeml
from pyrocko.io import stationxml as fdsn
from pyrocko.client import catalog
from pyrocko.automap import Map

from obspy.clients.fdsn.client import Client
from obspy import UTCDateTime
from obspy.core.event import Catalog
from obspy.core.stream import Stream
from obspy.core.event import Event
from obspy.core.event import Origin
from obspy.core.event import Magnitude
from obspy import read
from obspy import read_events
from obspy import read_inventory
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import pickle

import geopy.distance

#directory path
workdir='/Users/giaco/UNI/PhD_CODE/GIT/FLEGREI_moment_tensor'
datadir_raw=os.path.join(workdir,'DATA') 
datadir_sensitivity=os.path.join(workdir,'DATA_sensitivity')
meta_datadir=os.path.join(workdir,'META_DATA')

#read station.xml
stations_name=os.path.join(meta_datadir, 'stations_flegrei_INGV.xml')
stations=read_inventory(stations_name)
#print('stations selected:',stations)

#create dictionary with stations name and sensitivity values
st_sens={}
for network in stations:
    for stat in network.stations:
        st_sens[stat.code]={}
        for ch in stat.channels:
            if ch.code in st_sens[stat.code]:
                pass
            else:
                st_sens[stat.code][ch.code]=[]
            if ch.code=='HHE' or ch.code=='HHN' or ch.code=='HHZ':
                st_sens[stat.code][ch.code].append(ch.response.instrument_sensitivity.value)
#print(st_sens)

st_time={}
for network in stations:
    for stat in network.stations:
        st_time[stat.code]={}
        for ch in stat.channels:
            if ch.code in st_time[stat.code]:
                pass
            else:
                st_time[stat.code][ch.code]=[]
            if ch.code=='HHE' or ch.code=='HHN' or ch.code=='HHZ':
                st_time[stat.code][ch.code].append(ch.start_date)
                if ch.end_date!=None:
                    st_time[stat.code][ch.code].append(ch.end_date)
#print(st_time)

'''
#dictionary compr
d_apr = {'day_'+str(days) : [] for days in apr_days} 
# dizionario che contiene come keys i giorni da considerare per l'analisi
'''


'''
catname = os.path.join(catdir, 'catologue_flegrei_new_mag2_5.pf')

cat = model.load_events(catname)
print('Number of events:', len(cat))

client=Client('INGV')
stations_name=os.path.join(meta_datadir, 'stations_flegrei_INGV.xml')
stations=read_inventory(stations_name)                                 #read

print(stations)

datadir_3=os.path.join(workdir,'DATA_response')

for ev in cat:
    evID=ev.name

    #transform UTC time
    t = util.time_to_str(ev.time)

    print('origin UTC time event:',t)
    print('extimated magnitude:',ev.magnitude)

    event_start = UTCDateTime(t) - 20
    #print('event starts at:',event_start)

    event_end=UTCDateTime(t) +40
    #print('event ends at:',event_end)


    wave=Stream()
    for network in stations:
        for  station in network.stations:
            try:
                wave += client.get_waveforms(starttime=event_start,endtime=event_end,
                                    network=network.code,station=station.code,location='*', channel='HH?',
                                    attach_response=True)
            except:
                #print(station.code , 'station not recording')
                continue

    print('traces found:',len(wave.traces))

    wave.merge(fill_value=0)
    # trim over the [t1, t2] interval
    wave.trim(starttime=event_start, endtime=event_end, pad=True, fill_value=0)

    # remove trend
    wave.detrend("demean")
    
    #remove instrumental response
#    pre_filt = [0.02, 0.05, 25,30]       # for big eq
    pre_filt = [0.1, 0.2, 25,30]       # for small eq

    #remove instrumental response
    wave.remove_response(inventory=stations, output='DISP', pre_filt=pre_filt)

    waveletdir=os.path.join(datadir_3,evID)
    
    os.mkdir(waveletdir)

    wavelet_name= os.path.join(waveletdir,evID)  
    wave.write(wavelet_name +'.mseed',format='MSEED')
    print('saved!')
'''