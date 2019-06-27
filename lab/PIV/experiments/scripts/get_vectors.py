#--------------------------------------------------------------------------------------------------------
#
### Process all frame pairs (every 2 frames)
#
#--------------------------------------------------------------------------------------------------------
import f90nml
import numpy as np
import openpiv.tools
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import pandas as pd
from skimage import img_as_int
import sys
#--------------------------------------------------------------------------------------------------------

input_file = sys.argv[1] 
nml = f90nml.read(input_file) # namelist input parameters

# initial frame to process (int)
ini_frame = nml['frame_nml']['ini_frame']
# last frame to porcess (int)
end_frame = nml['frame_nml']['end_frame']
# pix/cm, get this scale from reference image
sca_factor = nml['frame_nml']['sca_factor'] 
# directory where frames are
frame_dir =  nml['frame_nml']['frame_dir']
# filename format of frame files
frame_format =  nml['frame_nml']['frame_format']
# directory to save velocity data files
vel_dir = nml['velocity_nml']['vel_dir']
 # filename format for velocity files
vel_format = nml['velocity_nml']['vel_format']

# Correlation processing parameters
window_size = nml['piv_nml']['window_size']
overlap = nml['piv_nml']['overlap']
dt = nml['piv_nml']['dt']
search_area_size = nml['piv_nml']['search_area']
thresh = nml['piv_nml']['thresh'] # signal to noise threshold

#--------------------------------------------------------------------------------------------------------

for ii in range(ini_frame,end_frame-2,2):

    fr_num_a = ('%04d' %ii)
    fr_num_b = ('%04d' %(ii+2))
    
    file_a = (frame_dir+frame_format %fr_num_a)
    file_b = (frame_dir+frame_format %fr_num_b)
    print(file_a)
    fr_a  = openpiv.tools.imread(file_a)
    fr_b  = openpiv.tools.imread(file_b)
    
    frame_a_int = img_as_int(1-fr_a)
    frame_b_int = img_as_int(1-fr_b)
    
    fra = frame_a_int[200:,500:1700].astype(np.int32) # crop images
    frb = frame_b_int[200:,500:1700].astype(np.int32)

    u, v, sig2noise = openpiv.process.extended_search_area_piv(fra, frb,
                                                               window_size, 
                                                               overlap=overlap, 
                                                               dt=dt, 
                                                               search_area_size=search_area_size, 
                                                               sig2noise_method='peak2peak')
    x, y = openpiv.process.get_coordinates(image_size=fra.shape, window_size=window_size, overlap=overlap)
    umask, vmask, mask = openpiv.validation.sig2noise_val(u, v, sig2noise, threshold = thresh)
    uval, vval, mask = openpiv.validation.global_val(umask, vmask, (-5000, 5000), (-5000, 5000))
    xsca, ysca, usca, vsca = openpiv.scaling.uniform(x, y, uval, vval, scaling_factor = sca_factor )
    openpiv.tools.save(xsca, ysca, usca, vsca, mask, (vel_dir + vel_format %fr_num_a))