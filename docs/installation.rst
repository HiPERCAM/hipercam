.. Installation guide, started 24/10/2017

.. include:: globals.rst

Getting the software
********************

Installation
============

The |hiper| pipeline is written in Python, specifically Python 3.X. The code
is stored on github and can be downloaded with::

    git clone https://github.com/HiPERCAM/hipercam.git

This will create a sub-directory called 'hipercam' wherever you run it. You
should move into that sub-directory where you will find a file README.rst
with further instructions and repeated verbatim below:

.. literalinclude:: ../README.rst
   :language: none

Updating
========

You should keep the directory 'hipercam' created by the first ``git clone``
command because then it makes updating simple. It becomes a matter of going
to the hipercam directory, pulling the changes from github, and 
re-installing::

  cd XXX/YYY/hipercam
  git pull
  python setup.py install --user

(the last line may vary according to how you install it).

Setting up
==========

In order to use the |hiper| pipeline for normal reduction, you just need the
commands in your path. Try typing ``hplot`` (a command) to see if they
are. Apart from this there are a few environment variables to know about.

Connecting to the server
------------------------

If using the pipeline in conjunction with the |hiper| instrument, and you
wish to use your own laptop rather than the dat reduction PC, then you need to
tell it where the server is which you do by setting an environment variable
called HIPERCAM_DEFAULT_URL. Typically I do this in my .cshrc file (C-shell,
use export with bash) with::

  setenv HIPERCAM_DEFAULT_URL 192.168.1.2:8007/

The value here is appropriate if you are connected to the |hiper| private
network via ethernet at the telescope. Note that there is no leading 'http://'
unlike the equivalent in the ULTRACAM pipeline because this address is used
with both http and web sockets and the appropriate prefix is tacked on by
whichever script is being used.  Also note the final '/' on the port number as
in '8007/', which must be there. When setting up for the first time, assuming
that the |hiper| file server is going, then the command
:func:`hipercam.scripts.hls` is a good one to start with to see that you are
set up correctly.

Matplotlib-based plots
----------------------

It seems hard to get set the fontsize of matplotlib plots consistently across
devices. Therefore the following environment variable should allow the fonts
to be scaled in size to your preference::

  setenv HIPERCAM_MPL_FSCALE 0.5

which would make them all a factor of 2 smaller.

