import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.interpolate import interp1d

def profiler_11jun(time_conduino):
    # Input
    filename = '/ocean/kramosmu/MultipleCanyons/lab/Conduino/scripts/profiler_tracking_11jun19.txt'
    ini = 19
    end = 85

    # Get timeseries 
    df = pd.read_csv(filename, delim_whitespace=True, header=1)
    time_p = np.array(df['t'][:])
    ydist = np.array(df['y'][:])
    top_pos = ydist[ini]
    top_pos2 = ydist[end]
    ini_time = time_p[ini]
    end_time = time_p[end]

    time = time_p[ini:end+1]-ini_time
    position = ydist[ini:end+1]-top_pos

    # Find interp function velocity of the profiler
    time_small = time[0:-1]
    vel = (position[1:]-position[:-1])/(time[1:]-time[:-1])
    func = interp1d(time_small, vel)
    vel_int = func(time_conduino)
    dist = np.zeros(np.shape(vel_int))
    
    for ii in range(1,len(time_conduino)):
        dist[ii] = (vel_int[ii]*(time_conduino[ii]-time_conduino[ii-1]))+dist[ii-1] 

    return(dist)

def profiler_19may(time_conduino):
    
    # Input
    filename = '/ocean/kramosmu/MultipleCanyons/lab/Conduino/Calibration/notebooks/mass_A.txt'
    ini = 28

    # Get timeseries 
    df = pd.read_csv(filename, delim_whitespace=True, header=1)
    time_p = np.array(df['t'][:])
    ydist = np.array(df['y'][:])
    top_pos = ydist[ini]
    top_pos2 = ydist[-1]
    ini_time = time_p[ini]
    end_time = time_p[-1]

    time = time_p[ini:]-ini_time
    position = ydist[ini:]-top_pos

    # Find interp function velocity of the profiler
    time_small = time[0:-1]
    vel = (position[1:]-position[:-1])/(time[1:]-time[:-1])
    func = interp1d(time_small, vel)
    vel_int = func(time_conduino)
    dist = np.zeros(np.shape(vel_int))
    
    for ii in range(1,len(time_conduino)):
        dist[ii] = (vel_int[ii]*(time_conduino[ii]-time_conduino[ii-1]))+dist[ii-1] 
    
    return(dist*100) # m to cm beacuse file is in meters

def profiler_04aug(time_conduino):
    
    # Input
    filename = '/ocean/kramosmu/MultipleCanyons/lab/Conduino/scripts/profiler_05aug19.csv'
    ini = 0

    # Get timeseries 
    df = pd.read_csv(filename, delim_whitespace=False, header=1)
    time_p = np.array(df['t'][:])
    ydist = np.array(df['y'][:])
    top_pos = ydist[ini]
    top_pos2 = ydist[-1]
    ini_time = time_p[ini]
    end_time = time_p[-1]

    time = time_p[ini:]-ini_time
    position = ydist[ini:]-top_pos

    # Find interp function velocity of the profiler
    time_small = time[0:-1]
    vel = (position[1:]-position[:-1])/(time[1:]-time[:-1])
    func = interp1d(time_small, vel)
    vel_int = func(time_conduino)
    dist = np.zeros(np.shape(vel_int))
    
    for ii in range(1,len(time_conduino)):
        dist[ii] = (vel_int[ii]*(time_conduino[ii]-time_conduino[ii-1]))+dist[ii-1] 
    
    return(dist) # file is in cm

