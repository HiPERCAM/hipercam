.. HiPERCAM windows format, started 16/02/2018

.. include:: ../globals.rst

.. |fig-windows| replace:: :numref:`fig-windows`

|hiper|'s windows format
************************

.. _fig-windows:

.. figure:: windows.png
   :scale: 100 %
   :alt: |hiper| windows
   :align: center

   Figure showing how the |hiper| windows are defined. The thin rectangles at
   the left and right of the imaging area represent the "pre-scan" (zero light)
   regions.

Although it is not usually necessary to know the details of |hiper|'s window
structure, it may help in somce cases. |fig-windows| illustrates the main 
points, which are:

  1. Each CCD is split into 4 quadrants labelled E, F, G, H. These have
     dimensions 1024 in X by 512 in Y, making the whole imaging area of each
     CCD 2048 by 1024 pixels.

  2. The outputs are located at the corners of the CCD (somewhat more in fact
     if one allows for an equivalent masked area). Parallel transfers are
     in the Y direction towards the respective output, while the serial
     register moves charges in the X direction towards the respective output.

  3. In windowed mode there can be either one or two windows per quadrant. The
     sizes of these windows are the same in all quadrants in sets of
     four. i.e. If you have one window per quadrant of dimension `WIN1 NX` by
     `WIN1 NY` in quadrant E, it will have the same dimensions in quadrants
     F, G and H, but you can have another window, dimensions `WIN2 NX` by
     `WIN2 NY`, that can be different from the first, but again will be
     repreated in all quadrants. The two windows must not overlap in the 
     Y direction.

  4. The location of a given set of windows is defined by parameter `YS` as in
     `WIN1 YS` indicating the Y-position of the first row to be read out
     relative to the quadrant's output. If it was actually the very first row
     of the quadrant that could be read, the `YS` would be set = 1. The value
     of `YS` is the same for all four windows of a given set. The X-locations
     of the windows are independently set by four parameters, one per
     quadrant, giving the position of the first column to be read out relative
     to the quadrant's output. These are called `XSE`, `XSF`, `XSG` and `XSH`
     with a `WIN1` or `WIN2` qualifier. e.g. `WIN1 XSE` says where the
     left-most pixels of the first window in quadrant E fall in terms of
     columns. This would = 1 if it lined up with the extreme left of the CCD.

  5. The windows can be be binned in both X and Y.

  6. There can be pre-scan pixels (see |fig-windows|) in the data. The
     pre-scans are 50 pixels (unbinned). If used, then the X-binning must
     divide into 50.
