.. Scripts docs created on 03/02/2018

HiPERCAM pipeline commands
**************************

The pipeline commands are all implemented as functions (and can be called as
such within scripts). This allows them to be included more naturally within
the overall package. The documentation below was generated automatically from
the in-code comments. See the documentation on hipercam.cline for how to
supply the arguments on the command line.

Each command appears therefore as a function, followed by a highlighted line
showing the arguments one can use on the command-line. Inputs in square
brackets such as `[source]` are hidden by default; those in brackets
e.g. `(plot)` may depend upon the value of earlier inputs. It is always safest
when first running a command simply to type its name and hit enter and let the
command itself prompt you for input. Many commands have hidden parameters that
can only be revealed by typing ``<command> PROMPT``. These are usually
parameters that rarely need changing, but you are sure sometimes to need to
alter them. A good starter, and a command with a large number of parameters,
both hidden and revealed, is ``rtplot``.


.. automodule:: hipercam.scripts
   :members:

