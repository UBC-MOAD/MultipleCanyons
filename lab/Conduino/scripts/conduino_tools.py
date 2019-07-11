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

