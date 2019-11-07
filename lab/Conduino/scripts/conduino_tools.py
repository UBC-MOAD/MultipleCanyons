import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt

def filter_freq(time, data, f):
    '''filter a particular frequency (f) from a timeseries (data). 
    As in harmonic analysis chapter 5.5   from Data Analysis Methods in 
    Physical Oceanography and amatlab snippet from Gonzalo
    
    Input:
    time, data are 1D arrays 
    f is a scalar
    
    Output:
    Yr: 1D array residual timeseries (Yr = data - Yh)
    Yh: 1D array with harmonic
    '''
    ine_freq = f
    x0 = np.ones_like(time)
    x1 = np.cos(2*np.pi*ine_freq*time)
    x2 = np.sin(2*np.pi*ine_freq*time)
    Yh = np.zeros_like(time)
    X = np.array([x0, x1, x2]) 

    D = np.matmul(X,np.transpose(X))
    Y = np.matmul(X,np.transpose(data-np.mean(data)))
    B = np.matmul(np.linalg.inv(D),Y), # regression coefficients.
    Yh = B[0][0] + B[0][1]*x1 + B[0][2]*x2
    Yr = data - Yh # Resiudual
    return Yr, Yh

def filter_timeseries(record, winlen=39):
    
    '''as in filter_timesies in salish sea tidetools.py without doodson option'''
    
    filtered = record.copy()
    record_length = record.shape[0]
    w = (winlen - 1) // 2
    weight = np.zeros(w, dtype=int)
    weight[:] = 1
    centerval = 1
    
    #Loop through record
    for i in range(record_length):
        
        # Adjust window length for end cases
        W = min(i, w, record_length-i-1)
        Weight = weight[:W]
        Weight = np.append(Weight[::-1], np.append(centerval, Weight))
        if sum(Weight) != 0:
            Weight = (Weight/sum(Weight))
        
        # Expand weight dims so it can operate on record window
        for dim in range(record.ndim - 1):
            Weight = Weight[:, np.newaxis]
        
        # Apply mean over window length
        if W > 0:
            filtered[i, ...] = np.sum(record[i-W:i+W+1, ...] * Weight, axis=0)
        else:
            filtered[i, ...] = record[i, ...]
    
    return filtered


def butter_lowpass(lowcut, fs, order=5):
    '''INPUT
       lowcut::float , frequency above which to filter signal (Hz)
       fs::float , sampling frequency (Hz)
       order::int, ortder of the filter
       OUTPUT
       b:
       a:
       '''    
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = butter(order, low, btype='low')
    return b, a


def butter_lowpass_filter(data, lowcut, fs, order=5):
    '''INPUT
       lowcut::float , frequency above which to filter signal (Hz)
       data::array , signal to filter
       fs:: sampling frequency in Hz
       order::int, ortder of the filter
       OUTPUT
       y::array, filtered signal of same size as data
       '''   
    b, a = butter_lowpass(lowcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def densP_29Mar19(reading):
    ''' Calibration from 29 mar 2019 in calibration_29mar19.ipynb. Returns density fitted using a 3rd deg polynomial.
    Input - reading::array
    Output - fitP::array of size [reading]'''
    fitP = 0.997378+(0.006040*reading)+(0.001648*reading**2)+(0.000105*reading**3)
    return(fitP)

def densANK1_29Mar19(read):
    '''Calibration from 29 March 2019 in calibration_29mar19.ipynb'''
    rho_1 = 0.997387+0.004844*read+0.000312*read**2+0.000204*read**3 
    return(rho_1)

def densANK2_29Mar19(read):
    '''Calibration  from 29 March 2019 in calibration_29mar19.ipynb '''
    rho_2 = 0.997311+0.006653*read+0.003429*read**2+-0.000041*read**3 # March 29, after knocking off
    return(rho_2)

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

def densP_12Jul19(read):
    ''' Calibration from 12 Jul 2019 in ANK_P20_probes_12jul19.ipynb
    Probes ANK and P20 calibrated independently - use if conduinos at 
    canyon heads and profiler were not connected simultaneously'''
    fitP = 0.998225+(0.002238*read)+(0.001219*(read**2))+(0.000168*(read**3))
    return(fitP)

def densANK1_12Jul19(read):
    ''' Calibration from 12 Jul 2019 in ANK_P20_probes_12jul19.ipynb
    Probes ANK and P20 calibrated independently - use if conduinos at 
    canyon heads and profiler were not connected simultaneously'''
    rho_1 = 0.998103+(0.004159*read)+(0.000582*(read**2))+(0.000166*(read**3)) 
    return(rho_1)

def densANK2_12Jul19(read):
    ''' Calibration from 12 Jul 2019 in ANK_P20_probes_12jul19.ipynb
    Probes ANK and P20 calibrated independently - use if conduinos at 
    canyon heads and profiler were not connected simultaneously'''
    rho_2 = 0.998142+(0.003846*read)+(0.000306*(read**2))+(0.000506*(read**3)) 
    return(rho_2)

def densP_09Aug19(read):
    ''' Calibration from 09 Aug 2019 in ANK_P20_probes_09aug19.ipynb
    '''
    fitP = 0.998003+(0.002483*read)+(0.005088*(read**2))-(0.000490*(read**3))
    return(fitP)

def densANK1_09Aug19(read):
    ''' Calibration from 09 Aug 2019 in ANK_P20_probes_09aug19.ipynb
    '''
    rho_1 = 0.997724+(0.003955*read)+(0.000772*(read**2))+(0.000136*(read**3)) 
    return(rho_1)

def errorANK1_09Aug19(read):
    '''Calculate the derivative drho(r)/dr of the calibration curve to calulate the 
    error of a density value or array as deltrho = drho(r)/dr delta r'''
    dummy = np.linspace(np.nanmin(read),np.nanmax(read),len(read))
    density = densANK1_09Aug19(dummy)
    del_rho = interp1d(dummy[1:],(density[1:]-density[:-1])/(dummy[1:]-dummy[:-1]), fill_value='extrapolate')
    return(del_rho(read))


def errorANK2_09Aug19(read):
    '''Calculate the derivative drho(r)/dr of the calibration curve to calulate 
    the error of a density value or array as deltrho = drho(r)/dr delta r'''
    dummy = np.linspace(np.nanmin(read),np.nanmax(read),len(read))
    density = densANK2_09Aug19(dummy)
    del_rho = interp1d(dummy[1:],(density[1:]-density[:-1])/(dummy[1:]-dummy[:-1]), fill_value='extrapolate')
    return(del_rho(read))

def densANK2_09Aug19(read):
    ''' Calibration from 09 Aug 2019 in ANK_P20_probes_09aug19.ipynb
    '''
    rho_2 = 0.997546+(0.005198*read)-(0.000409*(read**2))+(0.000470*(read**3))
    return(rho_2)

def densP_17Aug19(read):
    ''' Calibration from 17 Aug 2019 in ANK_P21_probes_17aug19.ipynb
    '''
    rho_P = 0.996552+(0.022331*read)+(-0.037047*(read**2))+(0.035005*(read**3))
    return(rho_P)

def densANK1_17Aug19(read):
    ''' Calibration from 17 Aug 2019 in ANK_P21_probes_017aug19.ipynb
    '''
    rho_1 = 0.997004+(0.008108*read)+(-0.002737*(read**2))+(0.000769*(read**3)) 
    return(rho_1)

def densANK2_17Aug19(read):
    ''' Calibration from 17 Aug 2019 in ANK_P21_probes_17aug19.ipynb
    '''
    rho_2 = 0.996888+(0.009153*read)+(-0.004260*(read**2))+(0.001388*(read**3)) 
    return(rho_2)

def errorANK1_17Aug19(read):
    '''Calculate the derivative drho(r)/dr of the calibration curve to calulate 
    the error of a density value or array as deltrho = drho(r)/dr delta r'''
    dummy = np.linspace(np.nanmin(read),np.nanmax(read),len(read))
    density = densANK1_17Aug19(dummy)
    del_rho = interp1d(dummy[1:],(density[1:]-density[:-1])/(dummy[1:]-dummy[:-1]), fill_value='extrapolate')
    return(del_rho(read))


def errorANK2_17Aug19(read):
    '''Calculate the derivative drho(r)/dr of the calibration curve to calulate 
    the error of a density value or array as deltrho = drho(r)/dr delta r'''
    dummy = np.linspace(np.nanmin(read),np.nanmax(read),len(read))
    density = densANK2_17Aug19(dummy)
    del_rho = interp1d(dummy[1:],(density[1:]-density[:-1])/(dummy[1:]-dummy[:-1]), fill_value='extrapolate')
    return(del_rho(read))

def densANK1_05Nov19(read):
    ''' Calibration from 05 Nov 2019 in ANK_probes_05nov19.ipynb
    '''
    rho_1 = 0.997970+(0.003838*read)+(0.003729*(read**2))+(-0.000199*(read**3)) 

    return(rho_1)

def densANK2_05Nov19(read):
    ''' Calibration from 05 Nov 2019 in ANK_probes_05nov19.ipynb
    '''
    rho_2 = 0.998538+(-0.002338*read)+(0.007345*(read**2))+(-0.000210*(read**3)) 

    return(rho_2)
