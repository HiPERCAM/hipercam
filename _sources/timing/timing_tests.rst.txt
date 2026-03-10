.. documentation of HiPERCAM timing tests

.. include:: ../globals.rst

.. |fig-0022-4| replace:: :numref:`fig-0022-4`
.. |fig-0022-5| replace:: :numref:`fig-0022-5`
.. |fig-0024-4| replace:: :numref:`fig-0024-4`
.. |fig-0024-5| replace:: :numref:`fig-0024-5`
.. |fig-0026-4| replace:: :numref:`fig-0026-4`
.. |fig-0026-5| replace:: :numref:`fig-0026-5`
.. |fig-0027-4| replace:: :numref:`fig-0027-4`
.. |fig-0027-5| replace:: :numref:`fig-0027-5`
.. |fig-0028-4| replace:: :numref:`fig-0028-4`
.. |fig-0028-5| replace:: :numref:`fig-0028-5`
.. |fig-0029-4| replace:: :numref:`fig-0029-4`
.. |fig-0029-5| replace:: :numref:`fig-0029-5`
.. |fig-0032-4| replace:: :numref:`fig-0032-4`
.. |fig-0032-5| replace:: :numref:`fig-0032-5`

Timing tests
************

To test the relations derived in the :doc:`timing` section, we rigged up an
LED pseudo-star on the focal plane mask. The LED turns on precisely and
near-instantaneously at the start of each UTC second (and off again about half
a second later). If one takes a series of short exposures (< 0.5 seconds),
then depending upon how they align with the square wave from the LED, one gets
either no counts, a steady maximum number of counts or something in between,
depending upon the degree of overlap between the exposure and the LED square
wave The data were extracted using |reduce|, and then plotted after folding on
the period of 1 second with an exact integer zero point. The data fold to give
either no counts, the full number of counts, and a linear ramp in between of
duration equal to the time spent accumulating photons. If the timing code is
correct, the mid-point of the ramp should coincide with zero, and any offset
gives an estimate of the precision to which times from |hiper| can be trusted.

WHT commissioning run, October 2017
===================================

This section groups a series of tests taken on 22 October 2017 during the
WHT commissioning run. The main issue as far as these timing tests was that
the voltage swings on the chips were a little too low leading to non-linear
behaviour. Some of the results need to be viewed as uncertain because of this.


run0022
-------

This was a no-clear windowed mode run from the first commissioning run on the
WHT in October 2017. |fig-0022-4| and |fig-0022-5| show the test data for
CCDs 4 and 5 (the others were either not operative or saturated).


.. _fig-0022-4:

.. figure:: 2017-10-22-0022-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0022 of 2017-10-22

   Timing test data for CCD 4, run0022 of 2017-10-22 shows a
   :math:`+165\,\mu\text{s}` offset. There are clear signs of non-linearity in
   the 'ramp' and the readout was known to exhibit significant non-linearity
   in this run that was subsequently fixed.

.. _fig-0022-5:

.. figure:: 2017-10-22-0022-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0022 of 2017-10-22

   Timing test data for CCD 5, run0022 of 2017-10-22 shows a
   :math:`-1\,\mu\text{s}` offset. There are clear signs of non-linearity in
   offset. The count levels are much lower in this arm, and non-linearity
   appears less than in CCD 4.

For this mode and run at least, the timing looks pretty good, albeit limited
by the unknown effects of non-linearity. The exposure time calculated for this
standard exposures in this mode (R + E, see :doc:`timing`) is 0.1614 secs,
compared to the measured values of 0.1659 and 0.1696 secs. This may again be
affected by non-linearity and needs checking.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 4     |   +165      |  0.1659    |   0.1614      |
+-------+-------------+------------+---------------+
| 5     |   -1        |  0.1696    |   0.1614      |
+-------+-------------+------------+---------------+


run0024
-------

This was a no-clear windowed mode run run with shorter exposures than run0022.
|fig-0024-4| and |fig-0024-5| show the data for
CCDs 4 and 5 (the others were either not operative or saturated).


.. _fig-0024-4:

.. figure:: 2017-10-22-0024-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0024 of 2017-10-22

   Timing test data for CCD 4, run0024 of 2017-10-22 shows a +8 microsecond
   offset. Counts and non-linearity appear lower than in |fig-0022-4|.

.. _fig-0024-5:

.. figure:: 2017-10-22-0024-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0024 of 2017-10-22

   Timing test data for CCD 5, run0024 of 2017-10-22 shows a -19 microsecond
   offset.

The measured exposure times of 0.00993 and 0.00980 secs compare with 
a calculated value of 0.00945 secs.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 4     |   +8        |  0.00993   |   0.00945     |
+-------+-------------+------------+---------------+
| 5     |   -19       |  0.00980   |   0.00945     |
+-------+-------------+------------+---------------+


run0026
-------

This was a no-clear windowed mode run run with shorter exposures than run0022.
|fig-0026-4| and |fig-0026-5| show the data for
CCDs 4 and 5 (the others were either not operative or saturated).


.. _fig-0026-4:

.. figure:: 2017-10-22-0026-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0026 of 2017-10-22

   Timing test data for CCD 4, run0026 of 2017-10-22.

.. _fig-0026-5:

.. figure:: 2017-10-22-0026-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0026 of 2017-10-22

   Timing test data for CCD 5, run0026 of 2017-10-22.

These show offset of -9 and -70 microsecs, and exposure lengths of
0.00500 and 0.00490 compared to a calculated 0.00458 seconds.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 4     |   +9        |  0.00500   |   0.00458     |
+-------+-------------+------------+---------------+
| 5     |   +80       |  0.00490   |   0.00458     |
+-------+-------------+------------+---------------+


run0027
-------

The key feature of this run is the use of the 'NSKIP' option. This makes the
exposures longer, and possibly induces more non-linear behaviour in the case
of CCD 4. The resulting offsets are quite large, but note that the exposure
are concomittantly long.

.. _fig-0027-4:

.. figure:: 2017-10-22-0027-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0027 of 2017-10-22

   Timing test data for CCD 4, run0027 of 2017-10-22.

.. _fig-0027-5:

.. figure:: 2017-10-22-0027-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0027 of 2017-10-22

   Timing test data for CCD 5, run0027 of 2017-10-22.

With NSKIP in operation, one can see a clear difference in the slopes
of the transitions in |fig-0027-4| versus |fig-0027-5|. The NSKIP for
CCD 5 in this case was 4 versus 3 for CCD 4.

Summary results:

+-------+-------------+------------+---------------+---------+
| | CCD | | Offset    | | Exposure | | Exposure    | | NSKIP |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] | |       |
+=======+=============+============+===============+=========+
| 4     |   -1736     |  0.06075   |   0.06406     |    3    |
+-------+-------------+------------+---------------+---------+
| 5     |   +226      |  0.08871   |   0.08801     |    4    |
+-------+-------------+------------+---------------+---------+


run0028
-------

Another examples woith NSKIP in use.

.. _fig-0028-4:

.. figure:: 2017-10-22-0028-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0028 of 2017-10-22

   Timing test data for CCD 4, run0028 of 2017-10-22.

.. _fig-0028-5:

.. figure:: 2017-10-22-0028-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0028 of 2017-10-22

   Timing test data for CCD 5, run0028 of 2017-10-22.

Summary results:

+-------+-------------+------------+---------------+---------+
| | CCD | | Offset    | | Exposure | | Exposure    | | NSKIP |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] | |       |
+=======+=============+============+===============+=========+
| 4     |   +742      |  0.04167   |   0.04010     |    2    |
+-------+-------------+------------+---------------+---------+
| 5     |   +57       |  0.04054   |   0.04010     |    2    |
+-------+-------------+------------+---------------+---------+


run0029
-------

Another examples with NSKIP in use. [Looks similar to run0027 as far
as I can tell.]

.. _fig-0029-4:

.. figure:: 2017-10-22-0029-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0029 of 2017-10-22

   Timing test data for CCD 4, run0029 of 2017-10-22.

.. _fig-0029-5:

.. figure:: 2017-10-22-0029-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0029 of 2017-10-22

   Timing test data for CCD 5, run0029 of 2017-10-22.

This is quite different from run0027, but rather few points have been
caught in transition, and one of them is quite discrepant. Definitely require
longer runs for these tests.

Summary results:

+-------+-------------+------------+---------------+---------+
| | CCD | | Offset    | | Exposure | | Exposure    | | NSKIP |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] | |       |
+=======+=============+============+===============+=========+
| 4     |   -429      |  0.06347   |   0.06406     |    3    |
+-------+-------------+------------+---------------+---------+
| 5     |   +1452     |  0.09158   |   0.08801     |    4    |
+-------+-------------+------------+---------------+---------+


run0032
-------

This is the first drift mode run. There were 3 junk windows in the masked area
at the start before real data emerged. There are 3 CCDs with valid data:

.. _fig-0032-3:

.. figure:: 2017-10-22-0032-3.png
   :scale: 100 %
   :alt: timing test data for CCD 3, run0032 of 2017-10-22

   Timing test data for CCD 3, run0032 of 2017-10-22.

.. _fig-0032-4:

.. figure:: 2017-10-22-0032-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0032 of 2017-10-22

   Timing test data for CCD 4, run0032 of 2017-10-22.

.. _fig-0032-5:

.. figure:: 2017-10-22-0032-5.png
   :scale: 100 %
   :alt: timing test data for CCD5, run0032 of 2017-10-22

   Timing test data for CCD 5, run0032 of 2017-10-22.

This is quite different from run0027, but rather few points have been
caught in transition, and one of them is quite discrepant. Definitely require
longer runs for these tests.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 3     |   -733      |  0.08123   |   0.08306     |
+-------+-------------+------------+---------------+
| 4     |   +293      |  0.08385   |   0.08306     |
+-------+-------------+------------+---------------+
| 5     |   -39       |  0.08359   |   0.08306     |
+-------+-------------+------------+---------------+


run0033
-------

This is the next and faster drift mode run. I have reduced all 4 available
CCDs which I think show again the effect of non-linearity (and possibly some
saturation in CCD 3).

.. _fig-0033-2:

.. figure:: 2017-10-22-0033-2.png
   :scale: 100 %
   :alt: timing test data for CCD 2, run0033 of 2017-10-22

   Timing test data for CCD 2, run0033 of 2017-10-22.

.. _fig-0033-3:

.. figure:: 2017-10-22-0033-3.png
   :scale: 100 %
   :alt: timing test data for CCD 3, run0033 of 2017-10-22

   Timing test data for CCD 3, run0033 of 2017-10-22.

.. _fig-0033-4:

.. figure:: 2017-10-22-0033-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0033 of 2017-10-22

   Timing test data for CCD 4, run0033 of 2017-10-22.

.. _fig-0033-5:

.. figure:: 2017-10-22-0033-5.png
   :scale: 100 %
   :alt: timing test data for CCD 5, run0033 of 2017-10-22

   Timing test data for CCD 5, run0033 of 2017-10-22.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 2     |   +314      |  0.03405   |   0.03322     |
+-------+-------------+------------+---------------+
| 3     |   -2047     |  0.02979   |   0.03322     |
+-------+-------------+------------+---------------+
| 4     |   +171      |  0.03395   |   0.03322     |
+-------+-------------+------------+---------------+
| 5     |   +35       |  0.03362   |   0.03322     |
+-------+-------------+------------+---------------+

run0034
-------

Another still faster drift mode run (CCD 3 was badly saturated).

.. _fig-0034-2:

.. figure:: 2017-10-22-0034-2.png
   :scale: 100 %
   :alt: timing test data for CCD 2, run0034 of 2017-10-22

   Timing test data for CCD 2, run0034 of 2017-10-22.

.. _fig-0034-4:

.. figure:: 2017-10-22-0034-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0034 of 2017-10-22

   Timing test data for CCD 4, run0034 of 2017-10-22.

.. _fig-0034-5:

.. figure:: 2017-10-22-0034-5.png
   :scale: 100 %
   :alt: timing test data for CCD 5, run0034 of 2017-10-22

   Timing test data for CCD 5, run0034 of 2017-10-22.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 2     |   +721      |  0.01713   |   0.01723     |
+-------+-------------+------------+---------------+
| 4     |   +225      |  0.01788   |   0.01723     |
+-------+-------------+------------+---------------+
| 5     |   +31       |  0.01773   |   0.01723     |
+-------+-------------+------------+---------------+


run0060
-------

This fast drift run appears to reveal an issue with the LED brightness
which does not seem to reach its maximum until around 3 msec after switch
on. This must have an effect upon the accuracy of these tests. I think
it particularly shows in the over-estimated exposure time since that is
calculated using the time to go from low to high, but the high value is
effectively too large. This could also add a few tenths of a millisec to the
offset.

.. _fig-0060-2:

.. figure:: 2017-10-22-0060-2.png
   :scale: 100 %
   :alt: timing test data for CCD 2, run0060 of 2017-10-22

   Timing test data for CCD 2, run0060 of 2017-10-22.

.. _fig-0060-4:

.. figure:: 2017-10-22-0060-4.png
   :scale: 100 %
   :alt: timing test data for CCD 4, run0060 of 2017-10-22

   Timing test data for CCD 4, run0060 of 2017-10-22.

.. _fig-0060-5:

.. figure:: 2017-10-22-0060-5.png
   :scale: 100 %
   :alt: timing test data for CCD 5, run0060 of 2017-10-22

   Timing test data for CCD 5, run0060 of 2017-10-22.

Summary results:

+-------+-------------+------------+---------------+
| | CCD | | Offset    | | Exposure | | Exposure    |
| | #   | | [|musec|] | |  [sec]   | | calc. [sec] |
+=======+=============+============+===============+
| 2     |   -35       |  0.00311   |   0.00263     |
+-------+-------------+------------+---------------+
| 4     |   +177      |  0.00334   |   0.00263     |
+-------+-------------+------------+---------------+
| 5     |   +328      |  0.00363   |   0.00263     |
+-------+-------------+------------+---------------+

