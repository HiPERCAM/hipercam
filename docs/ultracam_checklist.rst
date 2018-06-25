.. Links to checklists 21/06/2018

.. include:: globals.rst

ULTRACAM observing checklist
****************************

This checklist is designed with the NTT in mind; some only need doing at the
start of a run. If you are new to ULTRACAM or even experienced, please read
through this list from start to finish at least once, early in the run, if
only as a reminder. It contains a mix of statndard procedure and tips for
better observing practice.

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
     biases with the NTT shutter closed, all lights off and the focal plane
     mask in. Check your biases for problems by running |makebias| on each one.

  #. **Always** have the telescope spiralling when taking skyflats. Check that
     it is operating by looking at any stars in the field which should jump
     around.

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
