.. Scripts docs created on 03/02/2018

.. highlightlang:: rest

Available commands
******************

This page documents the pipeline commands which are all part of the
`hipercam.scripts` sub-module. Help on any given command can be obtained at
the terminal using e.g. ``pydoc hipercam.scripts`` or ``pydoc
hipercam.scripts.rtplot``. If you are new to the pipeline, it would be worth
reading the :ref:`command-calling` section on how to call the pipeline
commands and specify parameters.

Command parameters
------------------

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
-----------------------

The pipeline commands are distinct from standard unix commands in having a
'memory', which is implemented through storage of inputs in disk files for
each command. Parameters can be 'global' or 'local' depending upon whether
they are reset across multiple commands or just the command of interest. The
parameter memory along with the use of backslashes '\\' can save a huge amount
of typing making for efficient operation once you get up to speed.

