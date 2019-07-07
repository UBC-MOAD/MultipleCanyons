import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.interpolate import interp1d

def profiler_11jun(time_conduino):
    
    # Input
    filename = '/ocean/kramosmu/MultipleCanyons/lab/Conduino/profiler_tracking_11jun19.txt'
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

def densP_02May19(read):
    ''' Calibration from 02 May 2019 in ANK_P20_probes_02may19.ipynb'''
    fitP = 0.998176+(0.001940*read)+(0.001296*read**2)-(0.000073*read**3)
    return(fitP)

def densANK1_02May19(read):
    '''Calibration from 02 May 2019 in ANK_P20_probes_02may19.ipynb'''
    rho_1 = 0.998102+0.004567*read+0.000676*read**2+0.000056*read**3 
    return(rho_1)

def densANK2_02May19(read):
    '''Calibration from 02 May 2019 in ANK_P20_probes_02may19.ipynb '''
    rho_2 = 0.997984+0.004090*read+0.001643*read**2+0.000193*read**3 
    return(rho_2)

def densP_06Jul19(read):
    ''' Calibration from 06 Jul 2019 in ANK_P20_probes_06jul19.ipynb'''
    fitP = 0.998074+(0.001739*read)+(0.001236*read**2)+(0.000031*read**3)
    return(fitP)

def densANK1_06Jul19(read):
    '''Calibration from 06 Jul 2019 in ANK_P20_probes_06jul19.ipynb'''
    rho_1 = 0.998349+0.002087*read+0.001721*read**2+0.000050*read**3 
    return(rho_1)

def densANK2_06Jul19(read):
    '''Calibration from 06 Jul 2019 in ANK_P20_probes_06jul19.ipynb '''
    rho_2 = 0.997829+0.003691*read+0.000034*read**2+0.000560*read**3 
    return(rho_2)

