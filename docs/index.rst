.. HiPERCAM pipeline documentation master file, created by
   sphinx-quickstart on Tue Oct 24 13:03:06 2017.
   You can adapt this file completely to your liking, but
   it should at least contain the root `toctree` directive.

.. include:: globals.rst

The |hiper| pipeline manual
****************************

The |hiper| reduction software, known as the "pipeline" for short,
serves two purposes. First, it provides a means of displaying and reducing
|hiper|, ULTRACAM and ULTRASPEC data. Second,
it provides an "Application Programmers Interface" (API) to allow you
to access and manipulate the same sets of data. The first will be of
interest to anyone using |hiper|, ULTRACAM and ULTRASPEC; the second
to those who want to code their own scripts.  This manual covers both
aspects. It does not cover use of the |hiper| GUI or the finding chart
generator, hfinder, or the equivalent components for ULTRACAM and
ULTRASPEC. Details of these for |hiper| may be found on `github
<https://github.com/HiPERCAM>`_; and follows these links for much
more information on
`ULTRACAM <http://www.vikdhillon.staff.shef.ac.uk/ultracam/userman/userman.html>`_ and
`ULTRASPEC on the TNT <http://www.vikdhillon.staff.shef.ac.uk/ultraspec/ultraspec_tnt.html>`_.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   news
   observers
   commands
   api/api
   installation
   files
   api/bugs
   makingmovies
   phaseII

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
