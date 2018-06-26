.. News created on 26/06/2018

.. include:: globals.rst

|hiper| pipeline news
*********************

This page documents important changes in the pipeline, as a help to whether
you should update.

#. **26 June 2018, Version 0.16.XX**
   There have been very many recent changes in the pipeline as I have started
   to use it in earnest for reduction and have fixed a few problems
   encountered along the way. I would strongly advise anyone using the
   pipeline to update to the latest version. Above all I have improved the
   robustness of the aperture positioning making it less vulnerable to cosmic
   rays, removing an irritating overflow problem associated with Moffat
   profiles, and putting in several more checks to guard against positioning
   problems. Altogether these significantly reduce the chances that you lose
   apertures during a reduction.

   I have also implemented an experimental parameter ``fit_alpha`` which
   applies a fraction of the shift determined for an aperture relative to any
   references which is equivalent to smoothing the position differences.


