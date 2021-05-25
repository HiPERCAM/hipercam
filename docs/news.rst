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

#  :doc:`Updates since 1.0.0 <change_log_1.0.0_to_now>`. Most recent release: 1.0.0

#. **25 May 2021, Version 1.0.0**: fringe-correction added

     * fringe correction added to |nrtplot|, |grab|, |averun|, |reduce|. Associated
       new commands |makefringe| and |setfringe| added.
     * |nrtplot| improved responsiveness by updating plots while-u-wait
     * more synchronisation of default output names with inputs
     * |makemovie| for movie stills added
     * |joinup| to export joined up frames added
     * |rtplot| bug fix

   Here is the
   :doc:`full list of changes, 0.22.0 to 1.0.0 <change_log_0.22.0_1.0.0>`.

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

   Here is the
   :doc:`full list of changes, 0.21.0 to 0.22.0 <change_log_0.21.0_0.22.0>`.




