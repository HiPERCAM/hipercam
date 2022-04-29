"""Command line script to create a dark frame"""

import sys
import os
import time
import warnings
import signal

import numpy as np

from trm import cline
from trm.cline import Cline

import hipercam as hcam

__all__ = [
    "makedark",
]

############################################################################
#
# makedark -- combines frames of a run into one using clipped mean averaging,
# with bias subtraction. The exposure of the dark is corrected by whatever
# exposure applies to the bias.
#
#############################################################################


def makedark(args=None):
    """``makedark [source] (run first last [twait tmax] | flist) bias
    sigma output``

    Combines the frames of a single run into a single frame using clipped-mean
    averaging. This uses ``grab`` to get the frames and ``combine`` to combine
    them. It subtracts a bias and corrects the exposure time by the exposure
    time of the bias.

    Parameters:

       source  : string [hidden]
           Data source, four options:

             |   'hs' : HiPERCAM server
             |   'hl' : local HiPERCAM FITS file
             |   'us' : ULTRACAM server
             |   'ul' : local ULTRACAM .xml/.dat files
             |   'hf' : list of HiPERCAM hcm FITS-format files

           The standard start-off default for ``source`` can be set
           using the environment variable
           HIPERCAM_DEFAULT_SOURCE. e.g. in bash :code:`export
           HIPERCAM_DEFAULT_SOURCE="us"` would ensure it always
           started with the ULTRACAM server by default. If
           unspecified, it defaults to 'hl'.

       run : str [if source ends 's' or 'l']
           run name to access

       first : int [if source ends 's' or 'l']
           First frame to access

       last : int [if source ends 's' or 'l']
           Last frame to access, 0 for the lot

       twait : float [[if source ends 's' or 'l', hidden]
           time to wait between attempts to find a new exposure, seconds.

       tmax : float [[if source ends 's' or 'l', hidden]
           maximum time to wait between attempts to find a new exposure,
           seconds.

       flist : str [if source ends 'f']
           name of file list. Assumed that these are dias and dark corrected.

       bias : str
           Names of bias frame (made e.g. with |makebias|). This is
           so that the counts left in the dark frame are genuine dark
           counts which can be scaled by the ratio of exposure lengths
           during dark subtraction.

       sigma : float
           The value of 'sigma' to pass to the clipped mean combination in
           'combine'

       output : str
           name of final combined file

      .. Note::

         This routine writes the files returned by 'grab' to
         automatically generated files, typically in .hipercam/tmp, to
         avoid polluting the working directory. These are removed at
         the end, but may not be if you ctrl-C. You should check
         .hipercam/tmp for redundant files every so often

    """

    command, args = cline.script_args(args)

    # get inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("source", Cline.GLOBAL, Cline.HIDE)
        cl.register("run", Cline.GLOBAL, Cline.PROMPT)
        cl.register("first", Cline.LOCAL, Cline.PROMPT)
        cl.register("last", Cline.LOCAL, Cline.PROMPT)
        cl.register("twait", Cline.LOCAL, Cline.HIDE)
        cl.register("tmax", Cline.LOCAL, Cline.HIDE)
        cl.register("flist", Cline.LOCAL, Cline.PROMPT)
        cl.register("bias", Cline.LOCAL, Cline.PROMPT)
        cl.register("sigma", Cline.LOCAL, Cline.PROMPT)
        cl.register("output", Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        default_source = os.environ.get('HIPERCAM_DEFAULT_SOURCE','hl')
        source = cl.get_value(
            "source",
            "data source [hs, hl, us, ul, hf]",
            default_source,
            lvals=("hs", "hl", "us", "ul", "hf"),
        )

        # set a flag
        server_or_local = source.endswith("s") or source.endswith("l")

        if server_or_local:
            resource = cl.get_value("run", "run name", "run005")
            root = os.path.basename(resource)
            cl.set_default('output', cline.Fname(root, hcam.HCAM))
            first = cl.get_value("first", "first frame to grab", 1, 0)
            last = cl.get_value("last", "last frame to grab", first, 0)
            if last < first and last != 0:
                sys.stderr.write("last must be >= first or 0")
                sys.exit(1)
            twait = cl.get_value(
                "twait", "time to wait for a new frame [secs]",
                1.0, 0.0
            )
            tmax = cl.get_value(
                "tmax", "maximum time to wait for a new frame [secs]",
                10.0, 0.0
            )

        else:
            resource = cl.get_value(
                "flist", "file list", cline.Fname("files.lis", hcam.LIST)
            )
            first = 1

        bias = cl.get_value("bias", "bias to subtract", "bias")

        sigma = cl.get_value(
            "sigma", "number of RMS deviations to clip", 3.0, 1.0
        )
        output = cl.get_value("output", "output name", "bias")

    # Now the actual work.

    # We pass full argument lists to grab and combine because with None as the
    # command name, the default file mechanism is by-passed. 'prompt' is used
    # to expose the hidden parameters. Every argument needed is passed to
    # avoid any interactive prompting. Note that because of the Cline
    # mechanisms, in particular prompting of one parameter depending on
    # another, it can be tricky to get this right. The keyword 'list' can be
    # helpful in case of difficulty. i.e.  rather than just 'prompt' you would
    # use 'prompt', 'list'

    if server_or_local or bias is not None:

        print("\nCalling 'grab' ...")

        args = [None, "prompt", source, "yes", resource]

        if server_or_local:
            args += [str(first), str(last), str(twait), str(tmax)]

        args += [
            "no",
            "none" if bias is None else bias,
            "none",
            "none", "none", "f32",
        ]
        resource = hcam.scripts.grab(args)

    # 'resource' is a list of temporary files at this point

    with CleanUp(resource, server_or_local or bias is not None) as cleanup:

        # The above line to handle ctrl-c and temporaries
        if first == 1:
            # test readout mode if the first == 1 as, with non clear
            # modes, the first file is different from all others. A
            # warning is issued.
            with open(resource) as f:
                for line in f:
                    if not line.startswith('#') or not line.isspace():
                        first_frame = line.strip()
                        break

            mccd = hcam.MCCD.read(first_frame)
            instrument = mccd.head.get("INSTRUME", "UNKNOWN")
            if (
                    instrument == "ULTRACAM"
                    or instrument == "HIPERCAM"
                    or instrument == "ULTRASPEC"
            ):
                if "CLEAR" in mccd.head:
                    if not mccd.head["CLEAR"]:
                        warnings.warn(
                            "You should not include the first frame"
                            " of a run when making a dark from readout"
                            " modes which do not have clear enabled"
                            " since the first frame is different from"
                            " all others."
                        )
                else:
                    warnings.warn(
                        f"{instrument} has readout modes with both clears "
                        "enabled or not between exposures. When no clear "
                        "is enabled, the first frame is different from all "
                        "others and should normally not be included when "
                        "making a bias. This message is a temporary stop "
                        "gap until the nature of the readout mode has been "
                        "determined with respect to clears."
                    )

        print("\nCalling 'combine' ...")
        args = [
            None, "prompt", resource,
            "none", "none", "none",
            "c", str(sigma), "i",
            "yes", output,
        ]
        hcam.scripts.combine(args)

        # correct exposure time of dark frame by the exposure time of
        # the bias frame used
        dark = hcam.MCCD.read(cline.add_extension(output, hcam.HCAM))
        bias = hcam.MCCD.read(cline.add_extension(bias, hcam.HCAM))
        if "EXPTIME" in dark.head and "EXPTIME" in bias.head:
            dexpose = dark.head["EXPTIME"]
            bexpose = bias.head["EXPTIME"]
            dark.head["EXPTIME"] = dexpose - bexpose
            print(
                "Corrected dark exposure time from {:.4f} to {:.4f}".format(
                    dexpose, dexpose - bexpose
                )
            )
            dark.write(cline.add_extension(output, hcam.HCAM), True)
        else:
            warnings.warn(
                "Could not find exposure time (EXPTIME) in the dark and/or"
                " the bias hence could not correct it in the dark"
            )

        print("makedark finished")

class CleanUp:
    """
    Context manager to handle temporary files
    """
    def __init__(self, flist, temp):
        self.flist = flist
        self.temp = temp

    def _sigint_handler(self, signal_received, frame):
        print("\nmakedark aborted")
        sys.exit(1)

    def __enter__(self):
        signal.signal(signal.SIGINT, self._sigint_handler)

    def __exit__(self, type, value, traceback):
        if self.temp:
            with open(self.flist) as fp:
                for line in fp:
                    os.remove(line.strip())
            os.remove(self.flist)
            print('temporary files removed')
