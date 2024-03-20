import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from obspy import UTCDateTime
from obspy.core.event import Catalog
from obspy.core.event import Event
from obspy.core.event import Origin
from obspy.core.event import Magnitude
from obspy.core.event.base import Comment
from obspy.core.event import ResourceIdentifier
from obspy import read_events
from obspy import read_inventory
import pickle
from pyrocko import util

#read flegrei stations
workdir='/Users/giaco/UNI/PhD_CODE/GIT/FLEGREI_moment_tensor'
meta_data_dir= os.path.join(workdir,'META_DATA')

#%%% list of stations with GURALP CMG-40T-60S sensor
stations_xml_name=os.path.join(meta_data_dir,'stations_flegrei_INGV_original.xml')

inv_f=read_inventory(stations_xml_name)
#print(inv_f)

st_list=[]

for network in inv_f:
    for station in network:
        for channel in station.channels:
            if channel.sensor.description == 'GURALP CMG-40T-60S' and channel.data_logger.description == 'INGV GILDA':
                    trans_funct_type=channel.response.response_stages[0]._pz_transfer_function_type     #transfer function type
                    t1=util.time_to_str(channel.start_date)                   #start of recording time
                    try:
                        t2=util.time_to_str(channel.end_date)                     #end of recording time
                    except:
                        t2='None'
                    norm_factor=channel.response.response_stages[0].normalization_factor     #normalization factor: A0
                    g0=channel.response.response_stages[0].stage_gain                        #G0
                    try:
                        find_correct_stage=channel.response.response_stages[1].cf_transfer_function_type
                        if find_correct_stage == 'DIGITAL':
                            stage_gain=channel.response.response_stages[1].stage_gain       #stage gain
                        else:
                            stage_gain='!STAGE GAIN ERROR!'
                    except:
                        find_correct_stage=channel.response.response_stages[2].cf_transfer_function_type
                        if channel.response.response_stages[2].cf_transfer_function_type == 'DIGITAL':
                            stage_gain=channel.response.response_stages[2].stage_gain       #stage gain
                        else:
                            stage_gain='!STAGE GAIN ERROR!'
                    zeros=channel.response.response_stages[0].zeros                         #zeros
                    poles=channel.response.response_stages[0].poles                         #poles
                    sensitivity=channel.response.instrument_sensitivity.value               #sensitivity
                    st_list.append([[network.code +'_'+ station.code +'_'+ channel.code],[trans_funct_type],[t1],[t2],[sensitivity],[norm_factor],[g0],[stage_gain],[zeros],[poles]])

#st_list values:
#                    [0]network_station_channel [1]transfer_function_type
#                    [2]start_time  [3]end_time  
#                    [4]SENSITIVITY  [5]NORM_FACTOR  [6]G0  [7]STAGE_GAIN  [8]ZEROS  [9]POLES

#for ch in st_list:
#    print(ch[0],ch[2],ch[3])
#    print(ch[4],ch[5],ch[6],ch[7])
#    print(ch[8],ch[9])

save_list=False
if save_list:
    file_st_list=os.path.join(meta_data_dir,'stations_w_GURALP_40T_60S.txt')
    with open(file_st_list, "wb") as file:   
        pickle.dump(st_list,file)

#%%% correct stations values in xml file
stations_xml_name=os.path.join(meta_data_dir,'stations_flegrei_INGV.xml')

inv_f=read_inventory(stations_xml_name)
#print(inv_f)
        

for network in inv_f:
    for station in network:
        for channel in station.channels:
            if channel.sensor.description == 'GURALP CMG-40T-60S' and channel.data_logger.description == 'INGV GILDA':
                    trans_funct_type=channel.response.response_stages[0]._pz_transfer_function_type     #transfer function type
                    t1=util.time_to_str(channel.start_date)                   #start of recording time
                    try:
                        t2=util.time_to_str(channel.end_date)                     #end of recording time
                    except:
                        t2='None'
                    norm_factor=channel.response.response_stages[0].normalization_factor     #normalization factor: A0
                    g0=channel.response.response_stages[0].stage_gain                        #G0
                    try:
                        find_correct_stage=channel.response.response_stages[1].cf_transfer_function_type
                        if find_correct_stage == 'DIGITAL':
                            stage_gain=channel.response.response_stages[1].stage_gain       #stage gain
                        else:
                            stage_gain='!STAGE GAIN ERROR!'
                    except:
                        find_correct_stage=channel.response.response_stages[2].cf_transfer_function_type
                        if channel.response.response_stages[2].cf_transfer_function_type == 'DIGITAL':
                            stage_gain=channel.response.response_stages[2].stage_gain       #stage gain
                        else:
                            stage_gain='!STAGE GAIN ERROR!'
                    zeros=channel.response.response_stages[0].zeros                         #zeros
                    poles=channel.response.response_stages[0].poles                         #poles
                    sensitivity=channel.response.instrument_sensitivity.value               #sensitivity