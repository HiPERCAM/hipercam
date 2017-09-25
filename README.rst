hipercam
========

hipercam is the reduction package for the multi-CCD camera HiPERCAM built
through an ERC Advanced Grant awarded to Vik Dhillon. HiPERCAM is a multi-CCD
high-speed camera for astrophysical research.

Installation
============

hipercam is written in Python3; it does not support Python2.x It relies on the
following third-party packages. At the moment I suggest you get these
installed first; in the future I might make it more automatic::

  astropy    : widely used astronomical Python package with lots of useful
               stuff.

  matplotlib : plotting package

  numpy      : numerical data package

  trm.pgplot : my own Cythonised wrapper for PGPLOT which I wrote specifically
               for hipercam although it is entirely indepenedent of it. PGPLOT
               is needed for some of the faster plots as matplotlib is way too
               slow. PGPLOT itself must be installed for this to work. You can
               get trm.pgplot from my github site at
               https://github.com/trmrsh/trm-pgplot

Once you have all these, then its the usual installation:

  python setup.py install

or

  python setup.py install --prefix=top_level_install_directory

or

  python setup.py --user

depending on whether you have root priveledges and how you have set things up.


Tom Marsh




