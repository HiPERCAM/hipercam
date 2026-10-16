.. organising HiPERCAM reduction

.. include:: globals.rst

File organisation
*****************

This page contains advice on how to organise |hiper| reduction in terms of
directory structure and file names.

Night-by-night directories
==========================

Since |hiper| file names repeat from night-to-night (i.e. you can get
'run0023.fits' on multiple nights), reduction should always be organised
into directories with one per night. e.g. somewhere there should be a 
top-level directory with a series of directories like so:

.. code-block:: console

  2018-01-02/
  2018-01-03/
  2018-01-04/
  2018-02-11/
  2018-02-12/
  2018-02-13/

(or if you prefer, you could have underscores '_' rather than hyphens
'-'). This is how I organise the raw data in Warwick, and I recommend
the same layout for reduction. Thus when starting on a new night's-worth of
data, the first steps are typically as follows:

.. code-block:: console

  mkdir 2018-01-03
  cd 2018-01-03
  ln -s <path_to_raw_data>/2018-01-03 data

then within the sub-directory 2018-01-03, all the raw data will be found
inside 'data', a soft-link to the raw data directory.

File naming within a night directory
====================================

Typically a given night has several biases, a flat field and possibly some dark
runs. I always stick to the run name when making single calibration frames out
of such runs. Thus the raw |hiper| run 'run0003.fits' (or equivalently
run003.xml, run003.dat for ULTRACAM) will be combined in some way to make a
single |hiper| frame called 'run0003.hcm' (or 'run003.hcm'). However, then I make a
soft link which defines its fundamental purpose. e.g.:

.. code-block:: console

  ln -s run003.hcm bias_1x1_cdd.hcm

which says this is a bias frame taken with 1x1 binning with cdd readout speed
(ULTRACAM). If you set up 'ls' right, soft links appear distinct from ordinary
files and 'ls -l' should tell you where it points. Soft links always carry the
danger of being broken of course, so one must take care not to move or delete 
the target file. It helps also always to use relative path names (as long as
one does not change the directory structure), so I commonly make links to bias
frames taken on later nights e.g. within 2018-01-03 I might write:

.. code-block:: console

  ln -s ../2018-01-04/bias_arsco.hcm

where in 2018-01-04, bias_arsco.hcm is itself a link to perhaps 'run011.hcm',
a bias run taken specifically to match the format of a run on 'AR Sco' the
previous night. (NB without a link name specified in the previous line, the
link is simply created with the target name, i.e. 'bias_arsco.hcm' in this case.)

Always start bias frames with 'bias\_' then add enough information to make it
obvious what it can be used for. The binning is important, hence the '1x1'
above, but readout speeds must match as well, hence the 'cdd'. I adopt the
convention that I don't say when biases are full frame (no 'ff') because
that's the default bias, so the bias above is full frame.

For flat fields, the main issue is what filters were used, so I call them
names like 'flat_ugr.hcm' or 'flat_ugi.hcm'. Sometimes calibration frames are
affected by issues, but maybe you have no choice bu to use what is
available. Then a name like 'flat_ugr_caution.hcm' might be warranted.

Reduction file names
====================

Once all set with calibrations, then for a given run, I use the root
name ('run023' or 'run0034' or whatever) for all files associated with
that run. The pipeline defines a set of standard extensions that will
be added to any given file to allow this. Thus you will end up with
multiple files of the form 'run023.XXX', with 'XXX' being the file
extension. This makes it much easier to keep track of what a given
file is associated with. Many of the pipeline scripts are designed
with this in mind, e.g. |averun| will set the default output name
according to the run name. A full set of such files and their likely
meanings (including a few ancilliary files not associated with
specific runs) are:

===========  =========== =======================================================================
Name         Type        Meaning
===========  =========== =======================================================================
run023.hcm   Binary FITS Representative image of the run, probably from its start using |averun|
run023.ape   JSON text   File of photometric apertures
run023.red   ASCII text  Reduction driver file for |reduce|
run023.log   ASCII text  Output log file from a run of |reduce|
run023.fits  Binary FITS FITS version of a |reduce| log file made with |hlog2fits|
defect.dft   JSON text   File of CCD defects
fringe.frng  JSON text   File of peak/trough pairs for fringe measurement
files.lis    ASCII text  Column of file names (source='hf' option in some commands)
===========  =========== =======================================================================

If you find such a set, within a directory, then assuming the pipeline version
used to make the log file still applies, you should be able to re-reduce using
the aperture ('.ape') and reduction ('.red') files. A run of |setaper| with
the hcm ('.hcm') and ('.ape') files will show the apertures selected. If the
reduction fails, then you may want to re-generate the reduction driver file
using |genred|; see also |rupdate| which attempts to bring old reduction files
up to date.

