.. Observer's guide created on 25/10/2017

At the telescope
****************

Here it is assumed that you are starting to observe and want to look at your
data. It is written for someone starting afresh who just wants to get going.

The first command to try is 'rtplot'. Try typing in a terminal as follows::

  rtplot

You should see something like::

  rtplot
  run - run name [run0068]:

which is prompting you to supply a name for the run. If this is your first
time, the suggested value 'run0068', which remains from the last time the
command was run, is likely to be of little use. So you need to type the name
of your run, which might be 'run0002'. If you don't know it, either look at
the 'hdriver' observing GUI to see the run number, or ctrl-C out of 'rtplot' 
(safe to do so in any command) and type instead::

  hls

This returns a list of runs; you probably want the last one of the list, which
will be the most recent. Note that all HiPERCAM runs appear as 'run####'.




