.. News created on 26/06/2018

.. include:: globals.rst

|hiper| pipeline news
*********************

This page documents important changes in the pipeline, as a help to
whether you should update. I have revised how I list updates here and
from now on intend just to add new items when the first or second
number in the version changes, hence I have scrubbed old items (full
lists of past changes can be seen on github in any case). NB the
version of these docs may reflect a development version subsequent to
the most recent release and hence will appear to be one click on.

#. :doc:`Updates since 1.2.0 <change_log_1.2.0_to_now>`: patches released
   since the most recent minor version update.

#. **29 March 2022, Version 1.2.0**: many minor upgrades and fixes

     * |pbands| introduced to help make summary plots during observing
     * |exploss| to help judge the importance of readnoise on a run
     * |ncal| routine for noise calibration
     * |logsearch| database search routine working OK.

   :doc:`Full list of changes, 1.1.0 to 1.2.0 <change_log_1.1.0_1.2.0>`

#. **04 June 2021, Version 1.1.0**: over-haul of genred, hpackage, WCS in joinup

     * |genred| simplified but made more configurable by allowing one to define
       most parameters via a previously edited template reduce file
     * |hpackage| new command for bundling up standard reduce products
     * |joinup| use telescope pointing data for HiPERCAM+GTC to add a WCS.
     * improved handling of cleanup of temporary files
     * |nrtplot| and |joinup| both now include dark and fringe subtraction.
     * |joinup| can generate ds9 "region" files from |hiper| aperture files

   :doc:`Full list of changes, 1.0.0 to 1.1.0 <change_log_1.0.0_1.1.0>`.

#. **25 May 2021, Version 1.0.0**: fringe-correction added

     * fringe correction added to |nrtplot|, |grab|, |averun|, |reduce|. Associated
       new commands |makefringe| and |setfringe| added. Initial z-band fringe map
       for |hiper| created (for more see :ref:`fringing <fringing>`  and
       :ref:`fringing files <fringing_files>`)
     * |nrtplot| improved responsiveness by updating plots while-u-wait. This uses
       a hidden parameter `tupdate` that might require tuning a little according to
       the connection, although it may turn out not to matter.
     * more synchronisation of default output names with inputs (|makeflat|, |makebias|,
       |setaper|)
     * |makemovie| for movie stills added.
     * |joinup| to export joined up frames added.
     * |rtplot| bug fix that was causing problem with percentile scalings with small
       windows.

   :doc:`Full list of changes, 0.22.0 to 1.0.0 <change_log_0.22.0_1.0.0>`.

#. **05 May 2021, Version 0.22.0**: matplotlib-based real time image plotter
   introduced

     * |nrtplot|, matplotlib-based version of rtplot. Allows panning
       and zooming during the plots, choice of colour map, and
       profile fits from multiple CCDs.
     * setaper, hplot, setdefect also now support colour maps
     * averun output name now synchronised with input
     * new environment variable HIPERCAM_DEFAULT_SOURCE introduced
     * makebias output name synchronised with input
     * lots of minor improvements to logging scripts

   :doc:`Full list of changes, 0.21.0 to 0.22.0 <change_log_0.21.0_0.22.0>`.




