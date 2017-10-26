hipercam
========

hipercam is the reduction package for the multi-CCD camera HiPERCAM built
through an ERC Advanced Grant awarded to Vik Dhillon. HiPERCAM is a multi-CCD
high-speed camera for astrophysical research.

Installation
============

hipercam is written in Python3; it does not support Python2.x It relies on the
following third-party packages. At the moment I suggest you get these
installed first (you may well have several of them already); in the future I
might make it more automatic::

  astropy    : widely used astronomical Python package with lots of useful
               stuff.

  cython     : C-extensions for Python. Widely used package used to interface
               to C-libraries and to enable faster code when critical.

  matplotlib : plotting package

  numpy      : numerical data package. You need at least version 1.10
               providing numpy.stack

  trm.pgplot : my own Cythonised wrapper for PGPLOT which I wrote specifically
               for hipercam although it is entirely independent of it. PGPLOT
               is needed for some of the faster plots as matplotlib is way too
               slow. PGPLOT itself must be installed for this to work. You can
               get trm.pgplot from my github site at
               https://github.com/trmrsh/trm-pgplot
               It's not to be confused with "ppgplot" which is a more
               hand-crafted version with some differences in the calls.

  websocket  : for talking to the hipercam server. You want the Websocket
               client library which you can get from:
               https://github.com/websocket-client/websocket-client

Once you have all these, then its the usual installation:

  python setup.py install

or

  python setup.py install --prefix=top_level_install_directory

or

  python setup.py --user

depending on whether you have root priviledges and how you have set things up.
If you have multiple versions of python, then you need to specify python3


Tom Marsh




