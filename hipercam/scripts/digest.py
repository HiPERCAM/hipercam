import sys
import os
import shutil
import time
import glob
import re
import subprocess
import getpass
import argparse

import numpy as np
from astropy.time import Time, TimeDelta

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = [
    "digest",
]

#######################################################
#
# digest -- ingests hipercam data for archival purposes
#
#######################################################


def digest(args=None):
    description = """Ingests hipercam/ultra(cam|spec) data from the telescope for
    archiving purposes.

    This should be of little interest to most users. It does what the
    scripts checker, import_data, make_log_dirs and make_derived_dirs
    did for ULTRACAM data, but all in one go to reduce the need for
    thought.

    HiPERCAM run data from the telescope comes in a standarized
    form. The purpose of this script is run a few checks and set up a
    different standardised directory structure for it. It must be run
    inside a directory 'hipercam', 'ultracam' or 'ultraspec' which
    should contain a sub-directory 'raw_data' which contains run
    directories of the form 2017-10 which themselves contain
    night-by-night directories of the form YYYY_MM_DD. In these, it
    expects to find files in these of the form 'run1234.fits', (always
    4 digits), a file called MD5SUM_YYYY_MM_DD (i.e. matching the
    directory date), and a file called YYYY_MM_DD_log.dat. It will
    stop if any of these are not found and it is then up to the user
    to do something about it. It will also check that every run file
    appears in the MD5SUM file. Assuming these initial checks are
    past, it will then run 'md5sum -c' to check that the file md5sums
    match. 

    This program uses standard unix command-line switches. Run with
    '-h' to see help. It is not uncommon for there to be more log
    entries that there are runs if people pre-populate the log file
    but never take the runs. This causes an error that can be skipped
    with the '-i' option. 'digest' also now takes a more relaxed approach
    to any directories grouped under 'Others' where it will try to carry
    out checks but ignore problems where possible.

    """

    username = getpass.getuser()
    if username != "phsaap":
        print("digest is for archiving HiPERCAM runs and not " "meant for general use")
        return

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        dest="full",
        action="store_true",
        help="to a full (slow) search for newly imported data directories of the form YYYY_MM_DD",
    )
    parser.add_argument(
        "-i",
        dest="ignore",
        action="store_true",
        help="ignore runs not mentioned in the log but found in the directory",
    )
    args = parser.parse_args()
    full_search = args.full
    ignore = args.ignore

    RAW = "raw_data"
    DERIVED = "derived_data"
    LOGS = "logs"

    # check that we are in the right place
    subdirs = [sdir for sdir in os.listdir() if os.path.isdir(sdir)]

    basename = os.path.basename(os.getcwd())
    if (
        (basename != "hipercam" and basename != "ultracam" and basename != "ultraspec")
        or RAW not in subdirs
        or DERIVED not in subdirs
        or LOGS not in subdirs
    ):
        print(
            "** digest must be run in a directory called "
            "'hipercam', 'ultracam' or 'ultraspec' which "
            f"has sub-directories '{RAW}', '{DERIVED}' and '{LOGS}'",
            file=sys.stderr,
        )
        print("digest aborted", file=sys.stderr)
        return

    # regular expression to match run directory, e.g. 2018-10 [2018, October],
    # 2020-21 [season 2020 to 2021], 2021-P106 [ESO period 106]
    rdre = re.compile("^(Others|\d\d\d\d-(\d\d|P\d\d\d))$")

    # regular expression to match input night directory, e.g. 2018_10_30
    ndre = re.compile("^\d\d\d\d_\d\d_\d\d$")

    # regular expression to match hipercam raw data file
    rhre = re.compile("^run\d\d\d\d\.fits$")

    # regular expression to match ultra(cam|spec) raw data file
    rure = re.compile("^run\d\d\d\.dat$")

    # get the run directories
    rdirs = [
        os.path.join(RAW, rdir)
        for rdir in os.listdir(RAW)
        if os.path.isdir(os.path.join(RAW, rdir)) and rdre.match(rdir)
    ]
    rdirs.sort()

    if not len(rdirs):
        print("\n** Found no run directories", file=sys.stderr)
        print("digest aborted", file=sys.stderr)
        return

    print("\nChecking run-by-run ...\n")

    for rdir in rdirs:

        # test whether we might have done this already to save time
        # unless
        if not full_search:
            rd = os.path.basename(rdir)
            lrdir = os.path.join(LOGS, rd)
            drdir = os.path.join(DERIVED, rd)
            if os.path.exists(lrdir) or os.path.exists(drdir):
                continue

        print(f"Run directory = {rdir}")

        # strict checking
        strict = rdir.find('Others') == -1

        # get the imported night-by-night directories. These are
        # distinguished from the processed directories by having names
        # like 2018_10_31 rather than 2018-10-31
        ndirs = [
            os.path.join(rdir, ndir)
            for ndir in os.listdir(rdir)
            if os.path.isdir(os.path.join(rdir, ndir)) and ndre.match(ndir)
        ]
        ndirs.sort()
        nds = [os.path.basename(ndir) for ndir in ndirs]

        if len(ndirs):
            print(f"\nFound night directories: {nds}")
        else:
            print(f"\n** No night directories in {rdir}")
            continue

        for ndir in ndirs:

            # run elementary checks that are fast
            night = os.path.basename(ndir)

            log = os.path.join(ndir, f"{night}_log.dat")
            if os.path.isfile(log):
                print(f"... found the log file = {log}")

                # extract runs from the log file
                lruns = []
                with open(log) as fin:
                    for line in fin:
                        try:
                            name = line.split()[0]
                            lruns.append(name)
                        except IndexError:
                            pass

            else:
                lruns = None
                print(f"** no log file called {log} found", file=sys.stderr)
                if strict:
                    print("digest aborted", file=sys.stderr)
                    return

            md5 = os.path.join(ndir, f"MD5SUM_{night}")
            if os.path.isfile(md5):
                print(f"... found the md5sum file = {md5}")

                # extract runs from the MD5SUM file
                mruns = []
                with open(md5) as fin:
                    for line in fin:
                        hash, name = line.split()
                        if not name.endswith(".old"):
                            mruns.append(name[: name.rfind(".")])

            else:
                mruns = None
                print(
                    f"** no md5sum file called {md5} found", file=sys.stderr
                )
                if strict:
                    print("digest aborted", file=sys.stderr)
                    return

            # compile list of runs in the directory
            if basename.startswith("hiper"):
                # HiPERCAM
                runs = [
                    run[:-5]
                    for run in os.listdir(ndir)
                    if os.path.isfile(os.path.join(ndir, run)) and rhre.match(run)
                ]
            else:
                # ULTRACAM / ULTRASPEC
                runs = [
                    run[:-4]
                    for run in os.listdir(ndir)
                    if os.path.isfile(os.path.join(ndir, run))
                    and rure.match(run)
                    and os.path.isfile(os.path.join(ndir, run[:-4] + ".xml"))
                ]

            if strict and set(lruns) < set(runs):
                # complain if the log does not cover some of the runs
                print(
                    f"** {log} has missing runs cf the directory",
                    file=sys.stderr,
                )
                print(f"Missing runs are: {set(runs) - set(lruns)}")
                if ignore:
                    print("ignoring problem and continuing.")
                else:
                    print("digest aborted", file=sys.stderr)
                    return


            if strict and set(mruns) != set(runs):
                print(
                    "The runs in the md5sum file do not match "
                    "those in the directory",
                    file=sys.stderr,
                )
                print(f"Runs not in common are: {set(runs) ^ set(mruns)}")
                print("digest aborted", file=sys.stderr)

            print("... found all the runs listed in the md5sum file")

        # Now the md5sum checks
        print("\nNow running md5sum checks; these might take a while.\n")

        for ndir in ndirs:
            print(f"Night directory = {ndir}")
            night = os.path.basename(ndir)

            md5 = os.path.join(ndir, f"MD5SUM_{night}")
            if os.path.exists(md5):
                with open(md5) as fin:
                    for line in fin:
                        hash, run = line.split()
                        fname = os.path.join(ndir, run)

                        if os.path.exists(fname):

                            # run the md5sum
                            output = subprocess.check_output(["md5sum", fname]).decode()
                            hashc = output.split()[0]
                            if hashc == hash:
                                print(f"md5sum of {fname} is OK")
                            else:
                                print(
                                    f"** md5sum of {fname} is NOT OK!",
                                    file=sys.stderr,
                                )
                                if strict:
                                    print("digest aborted")
                                    return

        # Now change the data structure into my standard form:
        #
        # 1) Rename the night directories, replacing underscores in
        #    directory with '-' signs.
        # 2) Move each one up one level
        # 2) Create a soft-link to it in the night directory
        # 4) Copy the hand log of form yyyy_mm_dd_log.dat to yyyy-mm-dd.dat

        # create an equivalent run directory in the log directory
        rd = os.path.basename(rdir)
        lrdir = os.path.join(LOGS, rd)

        if not os.path.exists(lrdir):
            os.mkdir(lrdir)
            print(f"mkdir {lrdir}")

        # make a link to the telescope file
        tfile = os.path.join("..", "..", RAW, rd, "telescope")
        link = os.path.join(lrdir, "telescope")
        if not os.path.exists(link):
            os.symlink(tfile, link)
            print(f"ln -s {tfile} {link}")

        for ndir in ndirs:
            night = os.path.basename(ndir)
            nnight = night.replace("_", "-")
            nndir = os.path.join(RAW, nnight)
            lndir = os.path.join(os.path.dirname(ndir), nnight)

            # move the data directory up a level _ to -
            os.rename(ndir, nndir)
            print(f"\nmv {ndir} {nndir}")

            # link back into the run directory
            nfile = os.path.join("..", nnight)
            os.symlink(nfile, lndir)
            print(f"ln -s {nfile} {lndir}")

            # copy the hand-written log file to a new version without
            # underscores or 'log' in the name
            oldlog = os.path.join(nndir, f"{night}_log.dat")
            if os.path.exists(oldlog):
                newlog = os.path.join(nndir, f"{nnight}.dat")
                shutil.copyfile(oldlog, newlog)
                print(f"cp {oldlog} {newlog}")

            # make an equivalent night directory for log stuff
            logndir = os.path.join(LOGS, nnight)
            os.mkdir(logndir)
            os.chmod(logndir, 0o755)
            print(f"mkdir {logndir}")

            # create a link to the log night directory in the log run directory
            nfile = os.path.join("..", nnight)
            link = os.path.join(lrdir, nnight)
            os.symlink(nfile, link)
            print(f"ln -s {nfile} {link}")

            # create a link to the hand-written log in the log directory
            if os.path.exists(oldlog):
                lname = f"{nnight}.dat"
                lfile = os.path.join("..", "..", nndir, lname)
                link = os.path.join(logndir, lname)
                os.symlink(lfile, link)
                print(f"ln -s {lfile} {link}")

            # create a link to the data directory in the log directory
            dfile = os.path.join("..", "..", nndir)
            link = os.path.join(logndir, "data")
            os.symlink(dfile, link)
            print(f"ln -s {dfile} {link}")
