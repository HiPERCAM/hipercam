#
# This is an example "reduce file" which defines the operation of the reduce
# script. Basically it consists of a series of sections each of which contains
# a number of parameters. This file explains the meaning of these parameters.
#

# the 'apertures' section defines how the apertures are re-positioned from
# frame to frame. If not static, apertures are re-positioned through a
# combination of a search near a start location followed by a 2D fit. Several
# parameters below are associated with this process

[apertures]
file                = aper                 # file of software apertures for each CCD
reposition_mode     = reference_plus_tweak # static, individual, individual_plus_tweak, reference_plus_tweak
search_half_width   = 11                   # for initial search around start location, unbinned pixels
search_smooth_fwhm  = 4                    # smoothing FWHM, binned pixels

fit_method       = gaussian         # gaussian or moffat
fit_beta         = 4.               # Moffat exponent
fit_fwhm         = 6.               # FWHM, unbinned pixels
fit_fwhm_min     = 2.               # Minimum FWHM, unbinned pixels
fit_half_width   = 21               # for fit, unbinned pixels

# Calibration frames and constants
[calibration]
bias =        # Bias frame, blank to ignore
dark =        # Dark frame, blank to ignore
flat =        # Flat field frame, blank to ignore
readout = 3.  # RMS ADU. Float or string name of a file
gain = 1.     # Gain, electrons/ADU. Float or string name of a file
