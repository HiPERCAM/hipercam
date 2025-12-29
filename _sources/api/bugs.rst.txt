.. HiPERCAM bug reports, started 21/05/2018

.. include:: ../globals.rst

.. _reporting_problems:

Reporting problems
******************

If you encounter problems that you think are to do with the software,
then *please* raise an "issue" on `github
<https://github.com/HiPERCAM/hipercam>`_ where all code associated
with |hiper| is made public. The advantage of this over e-mail
is that all messages connected with a given problem are linked
together, and it makes it easy to see what remains to be done. Try to
make it one issue-per-problem.  If you e-mail me, I might well simply
respond with "open an issue on github".

The rules of the game are:

  #. Know which version you are running. Either type ``pydoc hipercam`` and
     scroll to the bottom of the page or type in a terminal::

        python -c 'import hipercam; print(hipercam.version())'

  #. Make sure you are up to date, i.e. compare your version with the version
     at the top of this page which will be a number like
     0.19.9.dev0+g52a7aac.d20200226 or 0.20.1 (the latter means its the
     version corresponding exactly to release 0.20.1; the former means you
     are using one developing towards 0.19.9).

  #. Try to narrow down the circumstances under which the problem arises; make
     sure it's repeatable. Specify how the problem occurs, cut-and-paste any
     error messages that are produced.

  #. Be prepared to tar or zip up a set of files (e.g. a reduce file, an
     aperture file, etc) that cause the issue so I can try to repeat it.



