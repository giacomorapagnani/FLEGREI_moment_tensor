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


workdir='/Users/giaco/UNI/PhD_CODE/GIT/FLEGREI_moment_tensor'

plotdir =  os.path.join(workdir,'PLOTS')
plotdir =  os.path.join(plotdir,'AMP_DIST')

datadir=os.path.join(workdir,'DATA_big_eq')                                         #CHANGE

###################################


for file in os.listdir(datadir):
    #select event
    name = os.fsdecode(file)

    if name.startswith('.'): 
        continue
    else:
        ev_dir=os.path.join(datadir,name)
        ev_name=os.path.join(ev_dir,name + '.mseed')

        #select wavelet (obspy)  
        w=read(ev_name)
        #print('number of traces in event:',len(w))
        w.plot(outfile= name + '.pdf' )

        '''
        #SAVE FIGURE SWITCH
        save_fig=True

        # Creazione della figura e dei subplot
        fig, axs = plt.subplots(1, 1, figsize=(17, 11), sharex=False)

        # Plot per il primo subplot
        plt.title(name)
        axs.scatter(num.array(distance1),
                        num.array(hhe),
                        label='HHE', s=20, color='green')
        axs.scatter(num.array(distance2),
                        num.array(hhn),
                        label='HHN', s=20, color='orange')
        axs.scatter(num.array(distance3),
                        num.array(hhz),
                        label='HHZ', s=20, color='blue')
        axs.set_xscale("log")
        axs.set_yscale("log")
        axs.set_ylabel('Amplitude')
        axs.grid(True)
        axs.set_xlabel('Distance [km]')
        axs.legend()

        for i, txt in enumerate(channel1):
            axs.annotate(txt, (distance1[i], hhe[i]),color='tab:green',size=7)

        for i, txt in enumerate(channel2):
            axs.annotate(txt, (distance2[i], hhn[i]),color='tab:orange',size=7)

        for i, txt in enumerate(channel3):
            axs.annotate(txt, (distance3[i], hhz[i]),color='tab:blue',size=7)

        if save_fig:

            figname = os.path.join(plotdir, name + '_amplitude_vs_distance.pdf')
            if os.path.isfile(figname):
                os.remove(figname)

            plt.savefig(figname)

            figname_svg = os.path.join(plotdir, name + '_amplitude_vs_distance.svg')
            if os.path.isfile(figname_svg):
                os.remove(figname_svg)

            plt.savefig(figname_svg)
            print('Figure',figname.split('/')[-1],'saved!')
        '''
#        plt.show()
        plt.close()
