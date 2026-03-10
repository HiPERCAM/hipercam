.. Installation guide, started 24/10/2017

.. include:: globals.rst

Getting the software
********************

|hiper| README file
===================

Here is the |hiper| README file which is written in the same markup language
as the rest of these pages.

.. include:: ../README.rst

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
that the |hiper| file server is going, then the command |hls| is a good one to start with to see that you are
set up correctly. Correspondingly, try |uls| if running the |hiper| pipeline
with the ULTRACAM server.

Matplotlib-based plots
----------------------

I have sometimes found matplotlib plots simply do not appear, e.g. the mean
level plot in |makebias|. If you encounter this problem, have a look at your
matplotlibrc file, probably located in .config/matplotlib in your home
directory. Mine say ``backend : Qt4Agg`` at the top, and you may need
something similar, depending on your matplotlib distribution.


Environment variables associated with |hiper|
---------------------------------------------

Here is a summary table of all environment variables associated with
|hiper| and their purpose. If you are simply reducing data at your home
institute, you probably won't need to set any of these.

=======================  ==============================================================
Environment variable     Purpose
=======================  ==============================================================
HIPERCAM_DEFAULT_SOURCE  Sets the default "source". Primarily useful for the instrument
                         control computer.
HIPERCAM_DEFAULT_URL     URL of the |hiper| FileServer on the rack PC to access data
HIPERCAM_ENV             Directory for storage of default parameter files
                         (~/.hipercam by default). Temporary files are stored in a
                         sub-drectory of this called "tmp". Although attempts are made
                         to delete such files, it is worth checking every now an again
                         that it does not contain too much junk.
ULTRACAM_DEFAULT_URL     URL of the ULTRACAM FileServer on the rack PC to access data
=======================  ==============================================================
