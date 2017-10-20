# This is an example "reduce file" which defines the operation of the reduce
# script. Basically it consists of a series of sections each of which contains
# a number of parameters. This file explains the meaning of these parameters.
# The idea is that these are to large extent unchanging and it would be annoying
# to be prompted every time for them.

[general]
version = 2017-10-20 # must match date in reduce itself.

ldevice  = 1/xs  # PGPLOT plot device for light curve plots
lwidth   = 0     # light curve plot width, inches, 0 to let program choose
lheight  = 0     # light curve plot height, inches

idevice  = 2/xs  # PGPLOT plot device for image plots [if implot True]
iwidth   = 0     # image curve plot width, inches, 0 to let program choose
iheight  = 0     # image curve plot height, inches

# the next section defines how the apertures are re-positioned from frame to
# frame. Apertures are re-positioned through a combination of a search near a
# start location followed by a 2D fit. Several parameters below are associated
# with this process. If there are reference apertures, they are located first
# to give a mean shift. The search on non-reference apertures can then be
# tightened.

[apertures]
aperfile   = aper  # file of software apertures for each CCD

search_half_width_ref  = 11   # for initial search around reference aperture, unbinned pixels
search_half_width_non  = 5    # for initial search around non-reference aperture, unbinned pixels
search_smooth_fwhm     = 4    # smoothing FWHM, binned pixels

fit_method     = moffat     # gaussian or moffat
fit_beta       = 4.         # Moffat exponent
fit_fwhm       = 6.         # FWHM, unbinned pixels
fit_fwhm_min   = 2.         # Minimum FWHM, unbinned pixels
fit_fwhm_fixed = no         # Slightly faster not to fit the FWHM.
fit_half_width = 21         # for fit, unbinned pixels
fit_sigma      = 4          # rejection threshold for fits
fit_height_min = 50         # minimum height to accept a fit

# The next lines define how the apertures will be re-sized and how the flux
# will be extracted from the aperture. There is one line per CCD which starts
# with "CCD label =" and then is followed by
#
# variable | fixed : how the aperture size will be determined. If variable it
# will be scaled relative to the FWHM, so there needs to be a FWHM.
#
# normal | optimal : how the flux will be extracted. 'normal' means a straight
# sum of sky subtracted flux over the aperture. This is often the best choice.
# 'optimum' tries to weight more towards regions of high flux, this is only
# possible if you have fit the profile earlier.
#
# Then follows a series of numbers in three triplets, each of which is a scale
# factor relative to the FWHM for the aperture radius if the 'variable' option
# was chosen, then a minimum and a maximum aperture radius in unbinned pixels.
# The three triples correspond to the innermost target aperture radius, the 
# inner sky radius and finally the outer sky radius. 

[extraction]
1 = variable normal 1.7 6.0 30.0 2.5 30.0 30.0 3.0 40.0 40.0
2 = variable normal 1.7 6.0 30.0 2.5 30.0 30.0 3.0 40.0 40.0
3 = variable normal 1.7 6.0 30.0 2.5 30.0 30.0 3.0 40.0 40.0
4 = variable normal 1.7 6.0 30.0 2.5 30.0 30.0 3.0 40.0 40.0
5 = variable normal 1.7 6.0 30.0 2.5 30.0 30.0 3.0 40.0 40.0

# next lines determine how the sky background level is calculated. Note
# you can only set error = variance if method = clipped. 'median' should
# usually be avoided as it can cause noticable steps in light curves. It's
# here as a comparator.
[sky]
method = clipped    # 'clipped' | 'median'
error  = variance   # 'variance' | 'photon': first uses actual variance of sky
thresh = 3          # threshold in terms of RMS for 'clipped'

# Calibration frames and constants
[calibration]
crop = yes    # Crop calibrations to match the data
bias =        # Bias frame, blank to ignore
dark =        # Dark frame, blank to ignore
flat =        # Flat field frame, blank to ignore
readout = 3.  # RMS ADU. Float or string name of a file
gain = 1.     # Gain, electrons/ADU. Float or string name of a file

# the light curve plot (includes transmission & seeing as well)

[lcplot]
xrange  = 0         # maximum range in X to plot (minutes), <= 0 for everything
extend_xrange = 10  # amount by which to extend xrange, minutes.

# light curve panel (must be present). Mostly obvious, then a series of lines,
# each starting 'plot' which specify one light curve to be plotted giving CCD,
# target, comparison ('!' if you don't want a comparison), an additive offset,
# a multiplicative scaling factor and then a colour for the data and a colour
# for the error bar There will always be a light curve plot, whereas later
# elements are optional, therefore the light curve panel is defined to have
# unit height and all others are scaled relative to this.

[light]
linear  = yes       # linear vertical scale (else magnitudes): 'yes' or 'no'
yrange_fixed = no   # keep a fixed vertical range or not: 'yes' or 'no'
y1 = 0              # initial lower y value
y2 = 0.2            # initial upper y value
extend_yrange = 0.1 # fraction of plot height to extend beyond points when rescaling

plot = 1 1 2 0 1 red red    # ccd, targ, comp, off, fac, dcol, ecol
plot = 1 3 2 0 1 red red    # ccd, targ, comp, off, fac, dcol, ecol
plot = 2 1 2 0 1 green red  # ccd, targ, comp, off, fac, dcol, ecol
plot = 2 3 2 0 1 green red  # ccd, targ, comp, off, fac, dcol, ecol

# configures the position plot. Can be commented out if you don't want one
# but make sure to comment it out completely, section name and all parameters
# you can have multiple plot lines

[position]
height  = 0.5    # height relative to light curve plot
x_fixed = no     # keep X-position vertical range fixed
x_min   = -5     # lower limit for X-position
x_max   = +5     # upper limit for X-position
y_fixed = no     # keep Y-position vertical range fixed
y_min   = -5     # lower limit for Y-position
y_max   = +5     # upper limit for Y-position
extend_yrange = 0.2  # Vertical extension fraction if limits exceeded
plot    = 1 1 green red

# configures the transmission plot. Can be commented out if you don't want one
# but make sure to comment it out completely, section name and all parameters
# you can have multiple plot lines

[transmission]
height = 0.5      # height relative to the light curve plot
ymax   = 110      # Maximum transmission to plot (>= 100 to slow replotting)
plot = 2 2 green red # CCD, target, data colour, error color

# configures the seeing plot. Can be commented out if you don't want one but
# make sure to comment it out completely, section name and all parameters you
# can have multiple plot lines. Don't choose linked targets as their FWHMs are
# not measured.

[seeing]
height = 0.5          # height relative to the light curve plot
ymax = 1.999          # Initial maximum seeing
scale = 0.3           # Arcsec per unbinned pixel
extend_yrange = 0.2   # Y extension fraction if out of range
plot = 2 2 green red  # CCD, target, data colour, error colour
