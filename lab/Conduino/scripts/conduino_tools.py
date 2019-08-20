import numpy as np

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

def densANK2_09Aug19(read):
    ''' Calibration from 09 Aug 2019 in ANK_P20_probes_19aug19.ipynb
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

