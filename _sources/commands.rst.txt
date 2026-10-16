.. command details doc created on 03/02/2018

.. include:: globals.rst

.. One-liners [hence ol] for each command that comes up repeatedly

.. |ol-add|      replace:: add two frames
.. |ol-averun|   replace:: average a series of frames in a run
.. |ol-cadd|     replace:: add a constant to a frame
.. |ol-calsearch| replace:: find matching calibration runs
.. |ol-cdiv|     replace:: divide a frame by a constant
.. |ol-cmul|     replace:: multiply a frame by a constant
.. |ol-combine|  replace:: combine a list of frames
.. |ol-csub|     replace:: subtract a constant from a frame
.. |ol-div|      replace:: divide one frame by another
.. |ol-exploss|  replace:: plot reduce log readout noise loss factor
.. |ol-fits2hcm| replace:: convert foreign FITs files to hcm format
.. |ol-flagcloud| replace:: flag cloudy and junk data in a reduce log
.. |ol-ftargets| replace:: automatically find targets using "sep"
.. |ol-genred|   replace:: create a reduce file
.. |ol-grab|     replace:: split frames out of a run
.. |ol-hfilter|  replace:: filter a HiPERCAM image
.. |ol-hinfo|    replace:: lists information on a HiPERCAM image
.. |ol-hist|     replace:: plot a histogram of a frame
.. |ol-hlog2col| replace:: produce ASCII column data from a reduce log
.. |ol-hlog2fits| replace:: convert a reduction log file to FITS
.. |ol-hls|      replace:: list the runs on the |hiper| server
.. |ol-hpackage| replace:: bundles up reduce data files
.. |ol-hplot|    replace:: plot a frame
.. |ol-joinup|   replace:: joins windows into a single image per CCD
.. |ol-ltimes|   replace:: list times of a run
.. |ol-ltrans|   replace:: computes transforms to align frames
.. |ol-makebias| replace:: combine a run to make a bias frame
.. |ol-makedark| replace:: combine a run to make a dark frame
.. |ol-makeflat| replace:: combine a set of frames into a flat
.. |ol-makefringe| replace:: combine a set frames into a fringe map
.. |ol-makemovie| replace:: makes stills for movies from a run
.. |ol-mstats|   replace:: list stats of multiple frames from a run
.. |ol-mul|      replace:: multiply two frames
.. |ol-ncal|     replace:: noise calibration
.. |ol-nrtplot|  replace:: plot frames as they come in [matplotlib]
.. |ol-pbands|   replace:: plot |reduce| log as a lightcurve per CCD
.. |ol-plog|     replace:: flexible plot of a |reduce| log
.. |ol-redanal|  replace:: analyse a reduction log file
.. |ol-reduce|   replace:: carry out photometric reduction
.. |ol-rupdate|  replace:: updates old reduce files
.. |ol-rtplot|   replace:: plot frames as they come in [pgplot]
.. |ol-setaper|  replace:: define the photometric apertures
.. |ol-setdefect| replace:: define a file of CCD defects
.. |ol-setfringe| replace:: define peak/trough pairs for fringe measurement
.. |ol-shiftadd| replace:: combine frames after aligning based on reduced positions
.. |ol-splice|   replace:: splice two frames together
.. |ol-stats|    replace:: report statistics of a frame
.. |ol-sub|      replace:: subtract two frames
.. |ol-uls|      replace:: list the runs on the ULTRACAM server

.. highlight:: rest

|hiper| commands
****************

This page documents the pipeline commands which are all part of
:mod:`hipercam.scripts`. Help on any given command can be obtained at the
terminal using e.g. ``pydoc hipercam.scripts`` or ``pydoc
hipercam.scripts.rtplot``. This is often the most useful way to get
information quickly. If you are new to the pipeline, it would be worth reading
the :ref:`command-calling` section on how to call the pipeline commands and
specify parameters.

.. contents:: Contents
   :local:

Main commands
=============

This section lists the main |hiper| commands [#f1]_ ; see also the bottom of the page
for all the documentation on the commands in one long list. Clicking on a
command name in the table below will take you to the relevant section of this long
list. The table also indicates contexts in which each command is particularly
useful.

.. table::
   :widths: 10 30 8 8 8 8 8

   +--------------+----------------+----------+----------+---------+-----------+------------+
   | Command      |Purpose         |Observing |Reduction |Plots &  |Arithematic|Information |
   |              |                |          |          |analysis |           |            |
   +==============+================+==========+==========+=========+===========+============+
   | |add|        | |ol-add|       |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |averun|     | |ol-averun|    | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |cadd|       | |ol-cadd|      |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |cdiv|       | |ol-cdiv|      |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |cmul|       | |ol-cmul|      |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |combine|    | |ol-combine|   | Yes      | Yes      |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |csub|       | |ol-csub|      |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |div|        | |ol-div|       |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |exploss|    | |ol-exploss|   | Yes      |          |         |   Yes     | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |fits2hcm|   | |ol-fits2hcm|  |          | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |flagcloud|  | |ol-flagcloud| |          | Yes      | Yes     |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |ftargets|   | |ol-ftargets|  | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |genred|     | |ol-genred|    | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |grab|       | |ol-grab|      | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hfilter|    | |ol-hfilter|   |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hinfo|      | |ol-hinfo|     |          |          |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hist|       | |ol-hist|      |          |          |  Yes    |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hlog2col|   | |ol-hlog2col|  |          | Yes      | Yes     |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hlog2fits|  | |ol-hlog2fits| |          | Yes      |  Yes    |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hls|        | |ol-hls|       | Yes      |          |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hpackage|   | |ol-hpackage|  | Yes      | Yes      |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |hplot|      | |ol-hplot|     | Yes      | Yes      |  Yes    |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |joinup|     | |ol-joinup|    | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |ltimes|     | |ol-ltimes|    |          |          |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |makebias|   | |ol-makebias|  | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |makedark|   | |ol-makedark|  |          | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |makeflat|   | |ol-makeflat|  |          | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |makefringe| | |ol-makefringe|| Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |makemovie|  | |ol-makemovie| |          |          | Yes     |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |mstats|     | |ol-mstats|    |          |          |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |mul|        | |ol-mul|       |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |ncal|       | |ol-ncal|      | Yes      | Yes      |  Yes    |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |nrtplot|    | |ol-nrtplot|   | Yes      | Yes      |  Yes    |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |pbands|     | |ol-pbands|    | Yes      |          |  Yes    |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |plog|       | |ol-plog|      | Yes      |          |  Yes    |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |redanal|    | |ol-redanal|   |          | Yes      |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |reduce|     | |ol-reduce|    | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |rtplot|     | |ol-rtplot|    | Yes      | Yes      |  Yes    |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |rupdate|    | |ol-rupdate|   | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |setaper|    | |ol-setaper|   | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |setdefect|  | |ol-setdefect| | Yes      |          |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |setfringe|  | |ol-setfringe| | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |shiftadd|   | |ol-shiftadd|  | Yes      | Yes      |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |splice|     | |ol-splice|    |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |stats|      | |ol-stats|     |          |          |         |           | Yes        |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |sub|        | |ol-sub|       |          |          |         |   Yes     |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+
   | |uls|        | |ol-uls|       | Yes      |          |         |           |            |
   +--------------+----------------+----------+----------+---------+-----------+------------+

.. _command-calling:

Parameter specification
=======================

The pipeline commands are distinct from standard unix commands in having a
'memory', which is implemented through storage of inputs in disk files for
each command, and also in prompting you if you don't specify a parameter on
the command line. Parameters can be 'global' or 'local' depending upon whether
they are reset across multiple commands or just the command of interest. The
parameter memory, along with the use of backslashes '\\\\' to accept default
values can save a huge amount of typing making for efficient operation once
you get up to speed.

Basic parameter input
---------------------

The command |rtplot| has more parameters than most others and is a
good one to start with (see also its replacement |nrtplot|). Suppose
then that we have a raw |hiper| file, :file:`run0076.fits`, that we
want to plot. If we type 'rtplot' and follow the prompts, the first
few lines might be::

  rtplot
  run - run name [run0064]: run0076
  first - first frame to plot [10]: 1

This shows that the last time |rtplot| was invoked, it was used to look
as :file:`run0064.fits` starting with frame 10. Note that the extension
'.fits' is not needed: |hiper| makes significant use of extensions to
differentiate between different forms of file, all of which could be
associated with the same run and therefore start with the 'run0064' or whatever.
Alternatively we could have typed::

  rtplot run0076
  first - first frame to plot [10]: 1

which would accomplish the same, setting the run name according to its
position in the command-line arguments. For this to work, you need to be
certain of the ordering. If you are not sure of the order, but you are sure
that the parameter is called ``run``, then::

  rtplot run=run0076
  first - first frame to plot [10]: 1

will do, and the first two parameters could be similarly specified::

  rtplot first=1 run=run0076

Note that by naming the parameters, the order becomes immaterial. Now,
assuming that the command was completed, if you run |rtplot| again you
might get::

  rtplot
  run - run name [run0076]:
  first - first frame to plot [1]:
  ccd - CCD(s) to plot [0 for all] [3]: 3 4
  nx - number of panels in X [2]:
  bias - bias frame ['none' to ignore] [none]:
  msub - subtract median from each window? [False]: True

You will see that the defaults for ``run`` and ``first`` have been updated by
the previous invocations of the command. This the parameter "memory" referred
to earlier.  Now suppose that we wish to repeat the command without changing
any parameter.  Then a simple::

  rtplot \\

will do. The double backslash '\\\\\\\\' indicates "take the default value for
all remaining parameters". Two are needed because the shell will strip one of
them. If you let at least one prompting line be output then a single backslash
will accomplish the same thing::

  rtplot
  run - run name [run0076]: \

If you wish to change just one parameter, say ``msub``, then much typing can
be saved with::

  rtplot msub=no \\

(NB: Boolean True/False parameters like ``msub`` can be reset False with either
'no', 'false' or even just 'n'.) Combined with up and down arrow keys, and
commands can be quickly repeated or modified.


Hidden parameters and special keywords
--------------------------------------

A look through the :ref:`command-definitions` will reveal that many commands
have "hidden" parameters that are not usually prompted for. The idea behind
this is to reduce the level of prompting, particularly for those parameters
that hardly ever change. The values of all parameters can be revealed through
the use of a special keyword ``list```::

  rtplot list \\
  rtplot
  source = hl
  device = 1/xs
  width = 0.0
  height = 0.0
  run = run0076
  first = 1
  twait = 10.0
  tmax = 20.0
  ccd = 3 4
  nx = 2
  pause = 0.0
  plotall = False
  bias = none
  msub = True
  iset = p
  plo = 10.0
  phi = 99.7
  xlo = 600.0
  xhi = 950.0
  ylo = 100.0
  yhi = 400.0

This reveals parameters ``source``, ``device``, ``width`` and ``height`` that
were never prompted for in the commands of the previous section. Their values
can be changed by giving another special keyword ``prompt``:

  rtplot prompt
  source - data source [hs, hl, us, ul, hf] [hl]:
  device - plot device [1/xs]:
  width - plot width (inches) [0.0]: 10
  height - plot height (inches) [0.0]: 8
  run - run name [run0076]:

etc. Usually once changed in this way, the hidden parameters will keep the
new values, although there are some which will return to a standard default
on the basis that it is almost always what is wanted and it would be dangerous
to default to a different value.

In addition to ``list`` and ``prompt``, there is a third special keyword
``nodefs`` which means do not read or write any default parameter values from
files. If ``nodefs`` is specified, then all parameters values need to be spelt
out; it's use is inside scripts in circumstances where multiple instances are
being used, or when there is interactive work going on too, to avoid causing
problems by over-writing the default files. Avoid using these keywords in any other
context. e.g. you would be asking for trouble if you named your files ``list.hcm``
or ``prompt.hcm``.


Strings with spaces
-------------------

If you need to specify a string with a space in it, either let a command
prompt you and type as normal, or quote it on the command line, e.g.::

  rtplot run="silly file name" \\

Setting a string to be blank
----------------------------

If you want to set a string to a blank, then a carriage return won't work
when prompted for it since this will just keep whatever default is
set. Instead type '' (two single quotes) which will be interpreted as
implying an empty string.

Parameter ranges
----------------

Many numerical parameters have mimimum and/or maximum values set. Typing
'?' at the prompt should tell you what the values are. Typing 'min' or 'max'
will set the parameter to have the appropriate value.

Tab completion
--------------

File name input is helped by tab completion: start typing, hit <tab>
and if the file exists, the name might be completed.

Global vs local parameters
--------------------------

All |hiper| pipeline parameters fall into one of two classes, either being
'local' to a command or 'global' to multiple commands. The ``run`` parameter
of |rtplot| for instance also appears in |grab| and if you change it in
|rtplot|, it will be changed in |grab|. This is very useful when running
a series of commands on the same file as the commands almost 'know' what you want,
saving much typing.

Trouble-shooting parameter input
--------------------------------

It very rarely happens that the file I/O needed to read the default values
can get confused. This can happen e.g. if the data type of a parameter has
changed. If very odd things seems to happen when you try to start a command,
then you might want to track down where the default files are located. Usually
this will be in `.hipercam` in your home directory, unless the environment
variable `HIPERCAM_ENV` has been defined to re-direct where the default files
are stored. Once you have located them it is always safe to delete one or more
or all of the default files (end with .def). The worst that happens is that
the commands have lost the default values.

Errors reported by scripts
==========================

If you use the pipeline for any significant time, you will at some point
get a rather forbidding-looking error traceback such as::

  reduce run=run013 \\
  Traceback (most recent call last):
    File "/storage/astro1/phsaap/software/python/bin/reduce", line 11, in <module>
      load_entry_point('hipercam==0.19.9.dev22+g758df9d.d20200303', 'console_scripts', 'reduce')()
    File "/storage/astro1/phsaap/software/python/lib/python3.6/site-packages/hipercam/scripts/reduce.py", line 405, in reduce
      with spooler.data_source(source, resource, first, full=False) as spool:
    File "/storage/astro1/phsaap/software/python/lib/python3.6/site-packages/hipercam/spooler.py", line 323, in data_source
      return UcamDiskSpool(resource, first)
    File "/storage/astro1/phsaap/software/python/lib/python3.6/site-packages/hipercam/spooler.py", line 108, in __init__
      self._iter = ucam.Rdata(run, first, False)
    File "/storage/astro1/phsaap/software/python/lib/python3.6/site-packages/hipercam/ucam.py", line 637, in __init__
      Rhead.__init__(self, run, server)
    File "/storage/astro1/phsaap/software/python/lib/python3.6/site-packages/hipercam/ucam.py", line 152, in __init__
      udom = xml.dom.minidom.parse(run + '.xml')
    File "/warwick/desktop/2018/software/MPI/GCC/7.3.0-2.30/OpenMPI/3.1.1/Python/3.6.6/lib/python3.6/xml/dom/minidom.py", line 1958, in parse
      return expatbuilder.parse(file)
    File "/warwick/desktop/2018/software/MPI/GCC/7.3.0-2.30/OpenMPI/3.1.1/Python/3.6.6/lib/python3.6/xml/dom/expatbuilder.py", line 910, in parse
      with open(file, 'rb') as fp:
  FileNotFoundError: [Errno 2] No such file or directory: 'run013.xml'

Yikes! Although it looks awful, it simply reflects
the chain of function calls that led to the problem, an extremely
useful diagnostic of problems in Python code. Such tracebacks look a
bit ugly, but almost certainly, in most cases, including this one,
they are caused by incorrect parameter inputs. The one to look at is
probably the last line or two, which reveals in this case that an expected
file 'run013.xml' was not found; a directory listing would confirm this.

Such errors follow from a standard Python approach of not trying to
add endless checking code, but to let the code tell you what happened
when errors are encountered. This has the merit of being very
informative (and de-clutters code), but it can make it hard to
distinguish between an essentially trivial issue, such as a missing
file, and a genuine problem with the code. If however you find you
can't get round a problem and the error reported does not look
innocent, then it might be a time to :ref:`report the
problem<reporting_problems>`. We will do our best to provide fast
solutions to critical issues.

.. _command-definitions:

Command definitions
===================

This section contains documentation auto-generated from the code. This
is the same as is returned from clicking command names in the lists at
the top of this page or from using ``pydoc`` in a terminal (e.g. try
``pydoc hipercam.scripts.reduce``). Each command appears as a function
(an implementation detail), followed by a highlighted line showing the
parameters one can use on the command-line. Inputs in square brackets
such as ``[source]`` are hidden by default; those in round brackets
e.g. ``(plot)`` may or may not be prompted depending upon earlier
inputs. It is always safest when first running a command simply to
type its name and hit enter and let the command itself prompt you for
input. Many commands have hidden parameters that can only be revealed
by typing e.g. ``rtplot prompt``. These are usually parameters that
rarely need changing, but you are sure sometimes to need to alter
them.  See the :ref:`command-calling` section for details on how to
specify command parameters.

In the one-line descriptions below, ``run`` refers to a complete run,
containing multiple images, stored in a .fits file. ``frame`` refers to a
single image from a run as might be extracted using |grab|. These have file
extension '.hcm' to distinguish them, although they are also FITS-format files.

.. autoapifunction:: hipercam.scripts.add
.. autoapifunction:: hipercam.scripts.averun
.. autoapifunction:: hipercam.scripts.cadd
.. autoapifunction:: hipercam.scripts.cdiv
.. autoapifunction:: hipercam.scripts.cmul
.. autoapifunction:: hipercam.scripts.combine
.. autoapifunction:: hipercam.scripts.csub
.. autoapifunction:: hipercam.scripts.div
.. autoapifunction:: hipercam.scripts.exploss
.. autoapifunction:: hipercam.scripts.fits2hcm
.. autoapifunction:: hipercam.scripts.flagcloud
.. autoapifunction:: hipercam.scripts.ftargets
.. autoapifunction:: hipercam.scripts.genred
.. autoapifunction:: hipercam.scripts.grab
.. autoapifunction:: hipercam.scripts.hfilter
.. autoapifunction:: hipercam.scripts.hinfo
.. autoapifunction:: hipercam.scripts.hist
.. autoapifunction:: hipercam.scripts.hlog2col
.. autoapifunction:: hipercam.scripts.hlog2fits
.. autoapifunction:: hipercam.scripts.hls
.. autoapifunction:: hipercam.scripts.hplot
.. autoapifunction:: hipercam.scripts.hpackage
.. autoapifunction:: hipercam.scripts.joinup
.. autoapifunction:: hipercam.scripts.ltimes
.. autoapifunction:: hipercam.scripts.ltrans
.. autoapifunction:: hipercam.scripts.makebias
.. autoapifunction:: hipercam.scripts.makedark
.. autoapifunction:: hipercam.scripts.makeflat
.. autoapifunction:: hipercam.scripts.makefringe
.. autoapifunction:: hipercam.scripts.makemovie
.. autoapifunction:: hipercam.scripts.mstats
.. autoapifunction:: hipercam.scripts.mul
.. autoapifunction:: hipercam.scripts.ncal
.. autoapifunction:: hipercam.scripts.nrtplot
.. autoapifunction:: hipercam.scripts.pbands
.. autoapifunction:: hipercam.scripts.plog
.. autoapifunction:: hipercam.scripts.redanal
.. autoapifunction:: hipercam.scripts.reduce
.. autoapifunction:: hipercam.scripts.rtplot
.. autoapifunction:: hipercam.scripts.rupdate
.. autoapifunction:: hipercam.scripts.setaper
.. autoapifunction:: hipercam.scripts.setdefect
.. autoapifunction:: hipercam.scripts.setfringe
.. autoapifunction:: hipercam.scripts.shiftadd
.. autoapifunction:: hipercam.scripts.splice
.. autoapifunction:: hipercam.scripts.stats
.. autoapifunction:: hipercam.scripts.sub
.. autoapifunction:: hipercam.scripts.uls

.. [#f1] Several other commands (``aligntool``, ``atanalysis``, ``atbytes``,
         ``calsearch``, ``harchive``, ``hlogger``, ``hmeta``, ``logsearch``,
         ``makefield``, ``makedata``, ``pfolder``, ``redplt``, ``tanalysis``, ``tbytes``)
         are not documented here as they are of specialist usage.  Information
         on these is however available via ``pydoc``.
