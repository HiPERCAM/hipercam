.. News created on 26/06/2018

.. include:: globals.rst

|hiper| pipeline news
*********************

This page documents important changes in the pipeline, as a help to whether
you should update.

#. **05 May 2021, Version 0.22.0**: matplotlib-based real time image plotter
   introduced

     * |nrtplot|, matplotlib-based version of rtplot. Allows panning
       and zooming during the plots, choice of colour map, and
       profile fits from multiple CCDs.
     * setaper, hplot, setdefect also now support colour maps
     * averun output name now synchronised with input

   Here is the
   :doc:`full list of changes, 0.21.2 to 0.22.0 <change_log_0.21.2_0.22.0>`.

#. **17 March 2021, Version 0.21.2**: bug fix for a problem with saving
   parameter default files.

     * new environment variable HIPERCAM_DEFAULT_SOURCE introduced
     * makebias output name synchronised with input
     * lots of minor improvments to logging scripts
     * filtid, new command, aimed at potential filter ID from sky flats

   Here is the
   :doc:`full list of changes, 0.21.1 to 0.21.2 <change_log_0.21.1_0.21.2>`.

#. **05 March 2021, Version 0.21.1**: again a large number of updates and
   fixes in this version. Most are small or to do with data logging, but
   the cumulative effect is significant and some bugs have been fixed.
   Tseries in particular has undergone a bit of an overhaul and has
   expanded functionality with an improved handling of bitmasks. Some
   of these may break existing code based upon Tseries; apologies. Main
   changes:

     * makedark bug fixed (was not subtracting bias)
     * new routine 'flagcloud' added to flag reduce logs interactively
     * reduce argument order and handling improved
     * a couple of new formats added to fits2hcm
     * Tseries: much better handling of bad and flagged data
     * Tseries bin option debugged
     * plog now has an option of specifying a plot title
     * should now be possible to display images with inverted axes
     * fixed bug preventing 'splice' from working
     * fixed bug for handling ULTRASPEC drift mode

   Here is the
   :doc:`full list of changes, 0.20.0 to 0.21.1 <change_log_0.20.0_0.21.1>`.

#. **10 July 2020, Version 0.20.0**: a large number of updates and
   fixes in this version. Apart from fixes to problems that would have
   been obvious as they were show-stoppers, the main changes are:

     * greatly improved target location for small windows
     * fixed bug with target location and parallel processing
     * added option to limit maximum beta for Moffat fits
     * added various specialist routines to do with timing
     * started on a routine |ftargets| which uses 'sep' to automate reduction

   There are many changes here, and I expect there to be bugs as a result,
   but some of the changes are significant enough that it is better to
   release now rather than wait longer. See also the
   :doc:`full list of changes, 0.19.8 to 0.20.0 <change_log_0.19.8_0.20.0>`

#. **26 February 2020, Version 0.19.8**: there have been many changes since
   the 0.16 versions. Apart from many bugfixes, the main ones are:

     * dark subtraction implemented
     * significant speed-ups implemented,
     * trimming of row & columns for ULTRACAM,
     * better display level setting through region definition,
     * ability to ignore frames with bad times in reduce,
     * rtplot bias level warnings added,
     * psf_reduce added for crowded fields.
     * rupdate script added to update old reduce files

#. **26 June 2018, Version 0.16.XX**: there have been very many recent changes
   in the pipeline as I have started to use it in earnest for reduction and
   have fixed a few problems encountered along the way. I would strongly
   advise anyone using the pipeline to update to the latest version. Above all
   I have improved the robustness of the aperture positioning making it less
   vulnerable to cosmic rays, removing an irritating overflow problem
   associated with Moffat profiles, and putting in several more checks to
   guard against positioning problems. Altogether these significantly reduce
   the chances that you lose apertures during a reduction.

   I have also implemented an experimental parameter ``fit_alpha`` which
   applies a fraction of the shift determined for an aperture relative to any
   references which is equivalent to smoothing the position differences.


