import sys
import os
import re
import time
import hipercam as hcam

__all__ = [
    "atbytes",
]

###########################################################################
#
# atbytes -- strips timing bytes out of all runs in a series of directories
#
###########################################################################


def atbytes(args=None):
    """``atbytes``

    Specialist script to strip out and save all timing bytes of all runs in a
    series of directories of YYYY[_-]MM[_-]DD form which should be
    sub-directories of the directory the script is run from. It will
    create (if necessary) a sub-directory of each of these night
    directories called 'tbytes' containing the timing bytes data, one
    for each run of the night directory. It can handle ULTRACAM,
    ULTRASPEC and HiPERCAM data. The script was created in preparation
    for making in-place modifications of the timing data after the
    discovery of some timing issues in March 2020.

    In the case of ULTRASPEC (or ULTRACAM, but really only applies to
    ULTRASPEC in practice) data, the script will always attempt to
    access any "old" version of the form "run023.dat.old" as these are
    files containing the original timing data potentially with null
    timestamps in the case of ULTRASPEC, as ultimately I want to
    develop a script that corrects both the ULTRASPEC null timestamps
    as well as the spurious but genuine timestamps that seem to be at
    the root of the March 2020 problem.

    atbytes takes no arguments. Just run it from a directory
    containing night directories.

    """

    # Specific formats for night directories and runs within them
    ndre = re.compile("^\d\d\d\d[_-]\d\d[_-]\d\d$")
    rnre = re.compile("^run\d\d\d\d\.fits|run\d\d\d\.xml$")

    # Find list of night directories
    ndirs = [
        dname for dname in os.listdir(".") if ndre.match(dname) and os.path.isdir(dname)
    ]
    ndirs.sort()

    for ndir in ndirs:

        print("\nStarted on", ndir)

        # create sub-directory where all timing info will be saved
        tbytes_dir = os.path.join(ndir, "tbytes")
        if not os.path.exists(tbytes_dir):
            os.mkdir(tbytes_dir)

        # read all the runs in the night directory
        fnames = [fname for fname in os.listdir(ndir) if rnre.match(fname)]
        fnames.sort()

        # change working directory to the tbytes one
        os.chdir(tbytes_dir)
        print("Changed working directory to", tbytes_dir)

        for fname in fnames:
            if fname.endswith(".fits"):
                source = "hl"
                run = os.path.join("..", fname[:-5])
                tfile = fname[:-5] + hcam.TBTS
                args = [None, "prompt", source, run]
            else:
                source = "ul"
                run = os.path.join("..", fname[:-4])
                tfile = fname[:-4] + hcam.TBTS
                args = [None, "prompt", source, "yes", run]

            if os.path.exists(tfile):
                print(tfile, "already exists; will not re-make")
            else:
                try:
                    hcam.scripts.tbytes(args)
                except hcam.ucam.PowerOnOffError:
                    print("ignoring", run, "which is a Power On or Off")
                except FileNotFoundError:
                    print("ignoring", run, "as no data were found")
                except (ValueError, hcam.HipercamError, KeyError) as err:
                    print("ignoring", run, "error =", err)

        print("Finished", ndir)

        # change up the working directory to the tbytes one
        os.chdir("../..")
        print("Moved working directory up two levels")
