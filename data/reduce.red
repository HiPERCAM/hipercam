# This is an example "reduce file" which defines the operation of the reduce
# script. Basically it consists of a series of sections each of which contains
# a number of parameters. This file explains the meaning of these parameters.
# The idea is that these are to large extent unchanging and it would be annoying
# to be prompted every time for them.

# the 'apertures' section defines how the apertures are re-positioned from
# frame to frame. If not static, apertures are re-positioned through a
# combination of a search near a start location followed by a 2D fit. Several
# parameters below are associated with this process

# If there are reference apertures, they are located first to give a mean
# shift. The search on non-reference apertures can then be tightened.

[apertures]
aperfile           = aper              # file of software apertures for each CCD
reposition         = yes               # yes / no to re-position the apertures
search_half_width_ref  = 11            # for initial search around reference aperture, unbinned pixels
search_half_width_non  = 5             # for initial search around non-reference aperture, unbinned pixels
search_smooth_fwhm = 4                 # smoothing FWHM, binned pixels

fit_method     = moffat     # gaussian or moffat
fit_beta       = 4.         # Moffat exponent
fit_fwhm       = 6.         # FWHM, unbinned pixels
fit_fwhm_min   = 2.         # Minimum FWHM, unbinned pixels
fit_fwhm_fixed = yes        # Slightly faster not to fit the FWHM.
fit_half_width = 21         # for fit, unbinned pixels
fit_sigma      = 4          # rejection threshold for fits
fit_height_min = 50         # minimum height to accept a fit 

# Calibration frames and constants
[calibration]
coerce = yes  # Try to force the calibrations to match the data
bias =        # Bias frame, blank to ignore
dark =        # Dark frame, blank to ignore
flat =        # Flat field frame, blank to ignore
readout = 3.  # RMS ADU. Float or string name of a file
gain = 1.     # Gain, electrons/ADU. Float or string name of a file

[general]
version = 2017-10-05 # must match date in reduce itself.
