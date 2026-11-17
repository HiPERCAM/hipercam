.. Links to checklists 21/06/2018

.. include:: globals.rst

ULTRACAM observing checklist
****************************

This checklist is designed with the NTT in mind; some only need doing at the
start of a run. If you are new to ULTRACAM or even experienced, please read
through this list from start to finish at least once, early in the run, if
only as a reminder. It contains a mix of statndard procedure and tips for
better observing practice. It's long because it is meant to be fairly
comprehensive. Ask me (TRM) if you don't follow any of the points.

  #. Make sure you have access to the phase II pages. Check the
     finding charts and verbal summaries for problems. Contact the PI if you
     find any. Watch out for finders made for the wrong telescope. They should
     say NTT on them if made with Stuart Littlefair's software.

  #. Download the NTT-compatible catalogue from the phase II. Get it copied
     across to where the TO can see it. 

  #. Check when sunset is so that you know when to be ready. The time for
     skyflats is ~10 to 30 minutes post sunset. Check that that
     the TO knows about the telescope spiral script.

  #. Check which filters are in place.

  #. Always take some biases in day time to check things are OK. Take science
     biases with the NTT shutter closed, all lights off, the focal plane mask
     in, with CLEAR enabled if possible, and NBLUE=1. When taking biases for
     specific targets and runs, the windows, readout speed and binning factors
     must all match. A good way to do this is simply to retrieve the setup, so
     it is always a good idea to save the configuration of every run. Remember
     however to set the exposure delay to zero and NBLUE=1 if you do
     this. Check your biases for problems by running |makebias| on each
     one. Take biases during daylight hours; don't waste valuable observing
     time unless the dome is closed because of weather.

  #. **Always** have the telescope spiralling when taking skyflats. Check that
     it is operating by looking at any stars in the field which should jump
     around. You can take flats in 'fbb' high-speed readout because you should
     have enough counts to rise above readout. However, don't be tempted to
     take very short exposures (below 2 seconds), because the fram transfer
     time will impose a slight slope on the results that gets worse at short
     exposures.

  #. Download the most up-to-date defect file (see :doc:`files`). This is of
     great value when acquiring targets to avoid the main bad spots.

  #. When acquiring, use the defect file, but don't rely on it alone. Look at
     your objects, especially the main target carefully in an averaged frame
     (see |averun|). Be prepared to move by a few pixels to steer clear of bad
     features and poor columns. **Do this for every target**

  #. Always check your target names via SIMBAD if possible. Do so via a
     browser rather than the udriver GUI "Verify" button as the latter has
     a tendency to cause the GUI to hang. If your target is new and not in
     SIMBAD then please use a name ending in "JHHMMSS.SS[+-]DDMMSS.S". These
     strictures are to help later archiving of the data. NB There is *no*
     other positional information in the ULTRACAM headers so correct object
     identification is essential.

  #. Make frequent comments on the weather, particularly clarity (or not) and
     seeing in the observing log. Also if a given run has an eclipse say, or
     something else that is notable about it, say so in the log. These sorts
     of remarks are valuable: you may be talking to your future self wondering
     "what on Earth?". Notable events include big focus shifts, satellite /
     meteor trails, big position jumps, unexpected objects in the field,
     asteroids, clouds, highly variable seeing and the like.

  #. Check the count levels in your target and comparisons in all CCDs. Are
     they likely to saturate? Might they exceed the "perppering" levels?
     Rememeber, these questions are highly seeing-dependent so remember that
     the seeing could improve. The peppering levels are well below saturation,
     especially in the green and blue CCDs. It's unclear whether peppering
     affects stars in the same way that it does flat field frames, but it is
     probably better to avoid exceeding these levels if possible.

  #. Once you start observing, you should try to get the reduction going asap.
     The main reason is to focus using the seeing plot. Make sure you are
     plotting the seeing measured in all CCDs and either use your main target
     or one nearby to it in case of CCD tilts. If using the hipercam pipeline
     (should do so), if you identify your target first with |setaper|, it will
     be assigned aperture label '1' and be the one that |genred| selects for
     the seeing plot by default. Your aim should be to get the seeing values
     of all CCDs overlapping each other. On the NTT, the focus tends to drift
     and you should monitor it all the time, and correct it at the start
     normally.

  #. Once you have seen to the focus, double-check your main target for poor
     columns or pixels. Move it if need be. If you have to make a large move,
     then it is often best to stop and re-start the run.

  #. Keep an eye on the target positions. I know of a run taken in 2017 which
     drifted in position by more than 50 pixels with no apparent attempt to
     correct the drift or comment in the log. This was bad enough to lose one
     of the comparisons altogether. If you have to move the targets, but don't
     want to break the run, try to move by relatively small amounts in
     possible a few steps as this will make the reduction simpler. The main
     requirement is not to move in a huge jump between consecutive frames.

  #. If you see different symbols appear in the light curve plot from
     |reduce|, the count levels in one or more of your targets may be too
     high. Check with |rtplot| which you can run at the same time as |reduce|.

  #. You should also check for overly low counts. Essentially it is a matter
     of the number of counts per pixel at the peak of each target compared to
     the gain*readout**2. You want the peak number of counts (sky plus target)
     to be significantly the larger to avoid being readout noise limited. For
     ULTRACAM this means >> 1.1*2.8**2 = 9 counts. I usually aim for at least
     50 to 100 counts at peak if possible.

  #. Make sure to write an entry *one line per target* in the standard ASCII
     log file for AutoLogger. Try to write stuff that is not otherwise easily
     seen. i.e. "Run on IP Peg" is not much extra use if the target name is
     "IP Peg", but "Bad cloud after 20 minutes" is useful.

  #. At the start of the night, add a comment to say what the filters are. If
     you change them, say e.g. "changed from r to i". This is a very helpful
     layer of redundancy on top of the filter names.

  #. If something goes badly wrong with a run, put "Junk!" somewhere in the
     comment, probably the start. 
