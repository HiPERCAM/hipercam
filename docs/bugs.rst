.. HiPERCAM bug reports, started 21/05/2018

.. include:: globals.rst

Reporting problems
******************

If you encounter problems that you think are to do with the software (not
unlikely), then please raise an "issue" on `github
<https://github.com/HiPERCAM/hipercam>`_. The great advantage of this over
e-mail is that all messages connected with a given problem are linked
together, and it makes it easy to see what remains to be done. Try to make it
one issue-per-problem.  If you e-mail me, I might well simply respond with
"open an issue on github".

Rules of the game are:

  #. Know which version you are running. The following line should tell you
     this::

     python -c 'import hipercam; print(hipercam.version())

  #. Make sure you are up to date (compare your version with the latest docs
     for example.

  #. Try to narrow down the circumstances under which the problem arises; make
     sure it's repeatable.

  #. Be prepared to tar up the set of files that cause the issue so I can try
     to repeat the issue.



