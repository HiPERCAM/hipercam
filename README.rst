hipercam
========

hipercam is the reduction package for the multi-CCD camera HiPERCAM built
through an ERC Advanced Grant awarded to Vik Dhillon. HiPERCAM is a multi-CCD
high-speed camera for astrophysical research.

Installation
============

hipercam is written in Python3; it does not support Python2.x It relies on the
following third-party packages. At the moment I suggest you get these
installed first (you will very likely have several of them already); in the
future I might try to make it more automatic:

  astropy :
         astronomical Python package with lots of useful
         stuff.

  cython :
         C-extensions for Python. Widely used package used to interface
         to C-libraries and to enable faster code when critical.

  matplotlib :
         standard plotting package for Python.

  numpy :
         Python's numerical data package. Highly likely you will have
         it, but note that you need at least version 1.12 to provide
         `numpy.flip`.

  trm.pgplot :
         my own Cythonised wrapper for PGPLOT which I wrote specifically
         for hipercam although it is entirely independent of it. PGPLOT
         is needed for some of the faster plots as matplotlib is way too
         slow. PGPLOT itself must be installed for this to work. Once
         you have PGPLOT, you can get trm.pgplot from my github site
         with::

             git clone https://github.com/trmrsh/trm-pgplot

         trm.pgplot is not to be confused with "ppgplot" which although
         very similar, is a more hand-crafted version with some
         differences in the calls.

  websocket :
         for talking to the hipercam server. You want the Websocket
         client library which you can get from::

             git clone https://github.com/websocket-client/websocket-client

  sep :
         astronomical source extractor for the aligntool script

  setuptools_scm :
         package to manage version numbers from the git repository


Once you have all these, then get the hipercam pipeline software itself using::

  git clone https://github.com/HiPERCAM/hipercam.git

This will create a directory (** reference from below). Change to that 
and then it's the usual Python installation procedure:

  python setup.py install

or

  python setup.py install --prefix=top_level_install_directory

or

  python setup.py --user

or

  pip install . --user

or

  pip install . --user --upgrade

depending on whether you have root privileges and how you have set things up.
If you have multiple versions of python, then you need to specify python3, and
the last line should use pip3.

If you have already installed but want to update, then go to the directory
created at step ** above and type ``git pull`` and then re-install.

For more information see:

  * `The documentation
    <http://deneb.astro.warwick.ac.uk/phsaap/hipercam/docs/html/>`_

  * `The |hiper| pipeline github repository <https://github.com/HiPERCAM/hipercam>`_

Tom Marsh
