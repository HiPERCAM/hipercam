.. documentation of HiPERCAM's photometry

.. include:: globals.rst

.. |fig-optimal| replace:: :numref:`fig-optimal`

Photometry in |hiper|
*********************

|hiper|'s |reduce| implements two forms of flux extraction, labelled
either 'normal' or 'optimal'. Normal extraction is essentially a matter of
adding up the flux above background in the target aperture. Optimal extraction
refers to a weighted extraction designed to yield the best signal-to-noise for
background-limited targets `(Naylor 1998)
<http://adsabs.harvard.edu/abs/1998MNRAS.296..339N>`_. This page discusses some
details of the implementations and the advantages and disadvantages of the
various options, and summarises some other aspects of how things work.

.. contents:: Contents
   :local:

.. _aperture_positioning:

Aperture positioning
====================

A key problem faced by the pipeline is how to locate the targets from frame to
frame given that the telescope moves, and conditions vary. The usual option
in the |reduce| configuration file is the ``[apertures]`` section is
``location = variable``. It's worth understanding how this operates. In the
usual case one defines one or more targets as ``reference`` apertures in
|setaper|. The aperture movement inside |reduce| then proceeds through a two
step process as follows:

 #. First search for each reference target in boxes of half width
    ``search_half_width`` around the last valid position of each one. This
    search is carried out by smoothing the image and taking the location of
    whatever local maximum exceeds a pre-defined threshold
    (``fit_height_min_ref``) and lies closest to the last-measured position.
    the maximum value as the starting position for a 2D profile fit. This
    method is fairly robust against even bright cosmic rays as long as they
    lie further from the expected position than the target.  If reference
    targets are chosen to be bright and isolated, one can carry out broad
    searches which allow for very poor guiding. Following the search, 2D
    profile fits are carried out and the mean x,y shift relative to the
    starting positions calculated. If this stage fails (e.g.  because of
    clouds), then the rest of the frame is skipped on the basis that if the
    reference targets cannot be located, then no others will be either.

 #. Next, the positions of non-reference, non-linked apertures are
    determined. This is done through 2D profile fits starting from the shift
    determined from the reference targets. An search is carried out to check
    that a sufficiently high maximum (``fit_height_min_nrf``) exists, but is
    only actually used to change position if there are no reference stars on
    the basis that it is more reliable to trust the reference stars than a
    search on a faint (potentially) target.  The idea is that the mean shift
    from the reference targets should provide a good start. An extra parameter
    ``fit_max_shift`` can be used to control how far the profile fits are
    allowed to wander from the initial positions, and another ``fit_diff`` is
    effective at weeding out discrepancies when multiple reference stars are
    used.

The combination of the options available in |setaper| and the |reduce|
configuration file are a powerful means to track objects over thousands of
exposures in a row. The main risk of failure comes when a combination of
clouds and cosmic rays cause a target to be too faint to register while a
cosmic ray does. Multiple reference stars and careful use of ``fit_max_shift``
and ``fit_diff`` can help in such cases. Once the aperture positions are
determined, |reduce| moves onto extracting the flux.

Sky background estimation
=========================

The sky background is estimated from a circular annulus centred on each
target. The user can define the inner and outer radii. Inside |setaper| it is
possible to mask out obvious stars from the sky annulus for any target. To
some extent, the two methods implemented (median and clipped mean) should
eliminate contaminating stars, but they are not perfect and if you can see an
obvious star, it is best to mask it. You should, if possible, try to set the
sky annulus radii such that there is a much larger area of sky than in the
target aperture, i.e. :math:`R^2_\text{out}-R^2_\text{in} \gg
R^2_\text{targ}`, in order to minimise the uncertainty due to the sky
estimate. The downside of this is the possibility of picking up more nearby
stars, starting to overlap the sky annulus with the edge of the window
within which the target falls and increasing the size so much that
quadratic variations in the sky background become important. Overlap of the
sky aperture with the edge of the window is not desirable as one always wants
the sky to be symmetrical around the object to eliminate any gradients in the
sky.

Normal extraction
=================

As indicated above, normal extraction is essentially a case of adding all the
pixels inside the inner circular aperture defining the target. It is not quite
as simple as this because the pixels are square-shaped (or rectangular for
non-equal binning factors) and not compatible with a circle, thus something
needs to be done about the pixels at the edge. If they are simply included or
excluded on the basis of how far their centres are from the centre of the
target aperture, one can suffer from pixellation noise that occurs when pixels
instantaneously appear or disappear owing to tiny change of target
position. To reduce this I use the commonly-adopted soft-edge weighting scheme
whereby a pixel which has its centre exactly on the edge of the circle gets
a weight of 0.5, which is linearly ramped from 0 to 1 as the pixel approaches
the aperture over a length scale comparable to its size.

An important element of normal extraction is the radius to extract the target
flux from. In an ideal world with steady seeing, this would have a single
value. However seeing is very rarely steady, and can exhibit very variable
behaviour in a night, or even within minutes. Thus a commonly used option is
the 'variable' aperture option (set in the reduce file). This scales the
radius of the extraction aperture to a multiple of the FWHM of the stellar
profiles of the given frame. A single mean value of FWHM is established
through profile fits to the targets and thus the same extraction radius is
used for all targets. Useful values for the scale factor range from 1.5 to
2.5. Stellar profiles can have extended wings, so one should expect not to
include all the target flux in the aperture. Instead the hope is that the same
fraction of flux is missed from both target and comparison, which then cancels
when performing relative photometry. The exact choice of scale is not easy, as
it represents a balance between increasing statistical noise for large radii
versus increasing systematic noise for small radii where the fraction of flux
lost from the apertures becomes large. Bright targets tend to push the choice
towards large radii, while faint targets are better with small radii. The key
is to try different values and compare the results side by side.

Sometimes one might have no suitable comparison star, or only ones of very
different colour to the target such that they are bright in some bands but
faint in others. If the weather is photometric, then extraction of just the
target with a large aperture, reminiscent of the old days of photometry with
photomultipliers through 10 or 15 arcsec diameter apertures, may be in order.
It is also possible to use the comparison star on one band to correct for
short timescale fluctuations in another where the comparison may be faint, but
one should then typically use a larger than normal radius because the profiles
of different bands are unlikely to be the same, and one also then needs to
allow for differential atmospheric extinction over longer timescales.

A final rare possibility is to use a target as its own comparison. This can be
useful if an object has very different behaviour at different wavelengths but
no suitable comparison. A good example of this is the star V471 |nbsp| Tau
which shows an eclipse of a white dwarf that is deep at short wavelengths but
shallow at long wavelengths.

Profile fitting
===============

In order to measure the FWHM to set the aperture size, |hiper| uses one of two
possible forms of stellar profile, namely a 2D gaussian or "Moffat" profile,
in each case symmetric. The gaussian profile is described by the following
relation

.. math::

   f(r) = h \exp \left(- \alpha r^2\right),

where :math:`r` is the distance from the centre of the gaussian, :math:`h` is
its height, and :math:`\alpha` is defined by the FWHM of the profile according
to

.. math::

   \alpha = \frac{4 \ln 2}{(\text{FWHM})^2}.

Moffat profiles (`Moffat 1969
<http://adsabs.harvard.edu/abs/1969A%26A.....3..455M>`_) usually provide much
better fits to stellar profiles than gaussians. Moffat profiles are described
by

.. math::

   f(r) = \frac{h}{(1 + \alpha r^2)^\beta},

where now :math:`\alpha` is a function of both the FWHM and :math:`\beta`:

.. math::

   \alpha = \frac{4(2^{1/\beta} - 1)}{(\text{FWHM})^2}.

Moffat profiles tend to gaussians as :math:`\beta \rightarrow \infty`, anything
a gaussian can fit, they can. Their advantage is their more extended
wings which often give a better fit to real stellar profiles. However, they
also have a downside in that some profiles seem to fall off more steeply than
gaussians and when fitted with Moffat profiles, there is a tendency for the
value of :math:`\beta` to climb to very large values. This probably depends
upon the telescope and instrument, so should be the subject of
experimentation. In order to be integrable, they should have :math:`\beta > 1`.
The profile fitting is implemented in :mod:`hipercam.fitting` on top of the
Levenberg-Marquardt fitting routine :func:`scipy.optimize.leastsq`. It is
designed to be able to sub-pixellate, i.e. the profile is evaluated at
multiple points within pixels to allow for cases where the seeing becomes
small compared to the pixels. This is controlled by the parameter `fit_ndiv`
in the reduce file. It's use does come with a speed penalty owing to the
larger number of function evaluations.


Optimal photometry
==================

For faint targets a significant improvement can be made by adding in the flux
in a weighted manner, as discussed by `Naylor (1998)
<http://adsabs.harvard.edu/abs/1998MNRAS.296..339N>`_. As Naylor discusses,
rather than the precisely "optimal weights" which are proportional to the
expected fraction of target flux in a pixel divided by its variance, it is
better in practice to use weights simply proportional to the expected fraction
alone, as determined with some sort of profile fit. This is because this
approach is robust to poor profile fits, which are not uncommon
given the vagaries of seeing and guiding. This contrasts with the optimal
extraction used in spectroscopy where the multiple pixels in the dispersion
direction allow the development of good models of the spatial profile.

An example of the improvement possible with optimal photometry in shown in
|fig-optimal|.

.. _fig-optimal:

.. figure:: optimal.png
   :scale: 100 %
   :alt: optimal versus normal photometry
   :align: center

   Normal (blue) versus optimal (green) photometry of the same target relative
   to a much brighter comparison. The optimal photometry has been moved up by
   0.01 for clarity.

In this case a faint target was chosen, and there is indeed an obvious
improvement in signal-to-noise. The aperture radius in this case was set to
1.8 times the FWHM of the profiles, and, while a smaller size might yield a
higher signal-to-noise for the normal photometry, one would not want to cut
down very much because of loss of flux.

The choice of weights is designed so that one does not need the profile to be
perfectly modelled, but it does require that both target and comparison have
the *same* profile. If this is not the case, e.g. the CCDs have significant
tilt leading to a FWHM that varies with position, then optimal photometry can
suffer from systematic errors that make things worse. As usual, such problems
manifest themselves most on bright targets, so one should not simply choose
'optimal' over 'normal' by default.
