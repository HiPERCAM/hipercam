hipercam
========

hipercam is the reduction package for the multi-CCD camera HiPERCAM
built through an ERC Advanced Grant awarded to Vik Dhillon. HiPERCAM
is a multi-CCD high-speed camera for astrophysical research. The
hipercam pipeline can also be used to reduce data from the ULTRACAM
and ULTRASPEC high-speed cameras, and other instruments too.

More information about using the hipercam package can be found in
`the documentation <https://vikdhillon.staff.shef.ac.uk/hipercam/pipeline/index.html>`_.

Docker Installation
===================

The easiest way to install hipercam is to create a Docker image. This
is simple and cross-platform, but the image is large (around 2.8GB). To 
install via this route, please download the
`dockerfile <https://raw.githubusercontent.com/HiPERCAM/hipercam/refs/heads/master/data/hipercam.dockerfile>`_

+++++
Linux
+++++

To build the docker image from the docker file, run::

  docker build -t hipercam:latest -f hipercam.dockerfile .

And, to run the image, use::

  docker run -it --rm -v <local-path-to-some-data>:/home/hiperuser/data hipercam:latest

This will mount `<local-path-to-some-data>` in the virtual machine at
`/home/hiperuser/data`. You can then run the pipeline as normal.

+++++
MacOS
+++++

Installation on Apple Silicon Macs is a little more complicated. The instructions below assume you have a modern
MacOS/Apple Silicon system with Docker Desktop installed. 

To build the docker image from the docker file, run::

  docker build -t hipercam:latest --platform linux/amd64 -f hipercam.dockerfile .

And, to run the image, use::
       
  xhost +
  docker run -it -e DISPLAY=host.docker.internal:0 --platform linux/amd64 --rm -v <local-path-to-some-data>:/home/hiperuser/data hipercam:latest

This will mount `<local-path-to-some-data>` in the virtual machine at
`/home/hiperuser/data`. You can then run the pipeline as normal.

Manual Installation
===================

hipercam is written in Python3; it does not support Python2.x It
relies on multiple third-party packages, as described in the next
section.  At minimum, you should ensure that Cython and trm.pgplot are
installed.  Once you have, get the hipercam pipeline software itself
using::

  git clone https://github.com/HiPERCAM/hipercam.git

This will create a directory; if you are updating then a simple "git
pull" from within this directory should suffice. Change to this
directory and install with::

  pip install . --user

If you have multiple versions of python, then you might need to
specify pip3. If you want a self-contained install, I suggest looking
into using "pipenv".

I use the Qt5Agg backend to ensure that the important command
|setaper| works. Others may work too of course. I have this set as the
default in my .config/matplotlib/matplotlibrc configuartion file with
the line `backend: Qt5Agg`

Third-Party Modules
===================

Apart from Cython and trm.pgplot, I hope that most of the extras will
get automatically installed if necessary by pip. So, if you have
Cython and PGPLOT ready, you might as well try ``pip install
. --user`` here and now.

If something seems amiss, here are details of the third-party packages
which you can either install via pip or by looking for the packages in
your O/S package manager. e.g. under fedora, Cython appears as
``python3-Cython``.

  astropy :
         astronomical Python package with lots of useful stuff.

  Cython :
         C-extensions for Python. Widely used package used to interface
         to C-libraries and to enable faster code when critical. It is
         needed at the setup stage so it might have to be installed first
         rather than relying on pip finding it, although I could be wrong.

  fitsio :
         Provides fairly direct access to FITS through the cfitsio library.
         Has some features that the pure python implementation in astropy
         lacks.

  keyring:
         Used in 'logsearch' to store some passwords safely.

  matplotlib :
         standard plotting package for Python. You will also need
         a Qt backend. In opensuse for instance the package
         python3-matplotib-qt5 enables the Qt5Agg backend.

  numba :
        Another package along with Cython to support faster numerics. Uses
        "just-in-time" compilation of selected routines.

  numpy :
         Python's numerical data package. Highly likely you will have
         it.

  pandas :
         Widely used table-handling module, based on numpy, and again
         quite likely installed already. It is used for some routines
         to do with generation of spreadsheets and logs in the
         pipeline (see also xlsxwriter below).

  requests :
         http request module. It may well be installed already.

  scipy :
         Scientific software for Python, closely linked to numpy.

  sep :
         astronomical source extractor based on Bertin's source extractor.

  reproject:
         The reproject package implements image reprojection methods 
         for astronomical images. This is an optional dependency,
         but you will need it for `shiftadd` to work.

  setuptools_scm :
         package to manage version numbers from the git repository

  trm.pgplot :
         Cython-ised wrapper for PGPLOT which I wrote specifically for
         hipercam although it is entirely independent of it. PGPLOT is
         needed for some of the faster plots as matplotlib is way too
         slow (at least as standardly used; I am hoping to move
         towards matplotlib-only). PGPLOT itself (F77/C library) must
         be installed for this to work. Once you have PGPLOT, you can
         get trm.pgplot from my github site with::

             git clone https://github.com/trmrsh/trm-pgplot

         trm.pgplot is not to be confused with "ppgplot" which,
         although very similar, is a more hand-crafted version with
         some differences in the calls. Once you have cloned it, you
         can enter trm-pgplot and install with pip, but make sure to
         set the enviroment variable PGPLOT_PNG inside setup.py to
         "true" or "false" first, according to whether you installed
         the PNG drivers with PGPLOT.

  trm.cline :
         handles command line parameters. Available on PyPi.

  trm.utils :
         generally useful routines used at a few places. Available
	 from PyPi.
	 
  websocket-client :
         for talking to the hipercam server.

  xlsxwriter :
         if you want to use the logging scripts hlogger, the object
         search script logsearch, or build log database tools that
         output xlsx files. Since these are unusual, the software is
         designed to build without insisting on this module.

Contributing
===================
The hipercam pipeline is made for its users, and we welcome contributions of many kinds.
If you notice a bug, or want to request new functionality, please 
`raise an issue <https://github.com/HiPERCAM/hipercam/issues>`_ or consider submitting
a `pull request <https://github.com/HiPERCAM/hipercam/pulls>`_.

Further Information
===================

For more information see:

  * `The documentation
    <https://vikdhillon.staff.shef.ac.uk/hipercam/pipeline/index.html>`_

  * `The HiPERCAM pipeline github repository <https://github.com/HiPERCAM/hipercam>`_

Tom Marsh
