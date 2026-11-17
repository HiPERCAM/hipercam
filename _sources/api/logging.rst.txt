.. Details of logging the data

.. include:: ../globals.rst

Data logging
************

The |hiper| pipeline includes several commands associated with the
creation and mantenance of data logs. They are not of general
interest, but this section gives a blow-by-blow account of how one
deals with new data. All the commands should be run from a directory
called "raw_data", which itself is a sub-directory of a directory
named after the instrument, i.e. hipercam, ultracam or ultraspec.
The sequence of steps is then as follows:

#. When some new data have been taken, say on 2017-04-03 to
   2017-04-05, define a name for the run (i.e. a bunch of nights),
   which could be "2017-04" for a run in April 2017. Create a
   sub-directory of this name inside raw_data. cd to it, and create a
   file called "telescope" containing the name of the telescope (WHT,
   NTT, TNT or VLT).

#. Next download the data saved by Paul Kerry's end_of_night_tasks
   into this directory. The nights should have dates of the form
   2017_04_03 etc.

#. cd up the raw_data directory. Run the command |harchive| with no
   arguments.  (If you add further nights later, you should run it
   with the switch '-f' or it will skip the run.) This checks for a
   set of expected files, changes the directory structure a bit and
   makes links. There is one special case run called "Others" for
   which the checks are lighter touch since they are not core data of
   the vikcam team.

#. Run |redplt| which looks for reduction log files and makes
   png plots of them.
   
#. Run |hlogger|. This works through all runs,
   extracting header and timing data (which will be written to a
   sub-directory 'meta' of each night directory. They also attempt to
   identify the target using the name in the header. This is often
   rather tricky as there are a fair few errors in the logs, both
   human and instrumental. For the vikcam team data, it is essential
   that all errors are fixed at this stage through a combination of
   editing header data and updating the various files TARGETS,
   AUTO_TARGETS, SKIP_TARGETS and FAILED_TARGETS. The switch '-n' is
   useful for quickly re-doing a single night, but once done a final
   complete run of |hlogger| will be needed to get all files correct.

#. If you kill |hlogger|, then delete whatever file it was
   in the midst of creating, i.e. delete 'meta/times' and
   'meta/posdata' in the night it had reached, otherwise you will end
   with partial files. You will normally need to kill using ctrl-Z
   followed by kill, since ctrl-C gets trapped (something I should
   probably change).
