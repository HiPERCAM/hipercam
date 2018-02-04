.. Scripts docs created on 03/02/2018

.. highlightlang:: rest

HiPERCAM commands
*****************

This page documents the pipeline commands which are all part of the
`hipercam.scripts` sub-module. Help on any given command can be obtained at
the terminal using e.g. ``pydoc hipercam.scripts`` or ``pydoc
hipercam.scripts.rtplot``. If you are new to the pipeline, it would be worth
reading the :ref:`command-calling` section on how to call the pipeline
commands and specify parameters.

.. contents::

.. _command-definitions:

Command definitions
===================

This section contains documentation auto-generated from the code.  The
pipeline commands are all implemented as functions (and can be called as such
within scripts). Each command therefore appears below as a function, followed
by a highlighted line showing the parameters one can use on the
command-line. Inputs in square brackets such as `[source]` are hidden by
default; those in round brackets e.g. `(plot)` may or may not be prompted
depending upon earlier inputs. It is always safest when first running a
command simply to type its name and hit enter and let the command itself
prompt you for input. Many commands have hidden parameters that can only be
revealed by typing e.g. ``rtplot prompt``. These are usually parameters that
rarely need changing, but you are sure sometimes to need to alter them.
See the :ref:`command-calling` section for details on how to specify command parameters.

.. automodule:: hipercam.scripts
   :members:

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

The command ``rtplot`` has more parameters than any other and is a good to
start with. Suppose we have a raw HiPERCAM file, :file:`run0076.fits`, that
we want to plot, then::

  rtplot
  run - run name [run0064]: run0076
  first - first frame to plot [1]:

would be the first few lines we might see if we entered ``rtplot`` and
followed the prompts. Alternatively::

  rtplot run0076
  first - first frame to plot [1]:

would accomplish the same, setting the run according to its position in
the command-line arguments. For this to work, you need to be certain of
the ordering. If you are not sure of the order, but you are sure that 
the parameter is called ``run``, then::

  rtplot run=run0076
  first - first frame to plot [1]:

will do, and the first two parameters could be similarly specified::

  rtplot first=1 run=run0076

Note that by naming the parameters, the order becomes immaterial. Now let's
go a bit further with ``rtplot``::

  rtplot
  run - run name [run0076]: run0076
  first - first frame to plot [1]: 
  ccd - CCD(s) to plot [0 for all] [3]: 3 4
  nx - number of panels in X [2]: 
  bias - bias frame ['none' to ignore] [none]: 
  msub - subtract median from each window? [False]: True

You will see that the default for ``run`` has been updated by the previous 
invocations of the command. This the parameter "memory" referred to earlier.
Now suppose that we wish to repeat the command without changing any parameter.
Then a simple::

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
the use of a special keyword ``list``::

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
can be changed by giving another special keyword ``prompt``::

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
problems by over-writing the default files.


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


