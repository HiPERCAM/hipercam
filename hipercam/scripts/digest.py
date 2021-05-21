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

###########################################################
#
# digest -- ingests hipercam etc data for archival purposes
#
###########################################################


def digest(args=None):
    description = """Ingests hipercam/ultra(cam|spec) data from the telescope for
    archiving purposes.

    This should not interest most users.

    |hipercam|, ULTRACAM and ULTRASPEC data come in standardised
    forms.  The purpose of this script is run a few checks and set up
    a different standardised directory structure for it. It must be
    run inside a directory 'hipercam/raw_data', 'ultracam/raw_data' or
    'ultraspec/raw_data' which should contain run directories of the
    form 2017-10, 2018-P104 or Others, and which themselves contain
    night-by-night directories of the form YYYY_MM_DD. In these, it
    expects to find files of the form 'run1234.fits' (|hipercam|) or
    'run123.dat' / 'run123.xml' (ULTRACAM, ULTRASPEC), along with a
    file called MD5SUM_YYYY_MM_DD (i.e. matching the directory date),
    and a file called YYYY_MM_DD_log.dat. It will stop if any of these
    are not found and it is then up to the user to do something about
    it. It will also check that every run file appears in the MD5SUM
    file. Assuming these initial checks are past, it will then run
    'md5sum -c' to check that the file md5sums match, and finally
    re-structure into the standard form for the Warwick archive.

    This program uses standard unix command-line switches. Run with
    '-h' to see help. It is not uncommon for there to be more log
    entries that there are runs if people pre-populate the log file
    but never take the runs. This causes an error that can be skipped
    with the '-i' option. 'digest' also now takes a more relaxed approach
    to any directories grouped under 'Others' where it will try to carry
    out checks but ignore problems where possible because it is common
    for other observers not to conform to our standards.

    .. Note::

       See also: |hlogger|, |ulogger|, |hmeta| and |redplt|

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
    parser.add_argument(
        "-s",
        dest="skip",
        action="store_true",
        help="skip MD5SUM checks",
    )
    args = parser.parse_args()
    full_search = args.full
    ignore = args.ignore
    skip_md5 = args.skip


    cwd = os.getcwd()
    if os.path.basename(cwd) != "raw_data":
        print("** digest must be run in a directory called 'raw_data'")
        print("digest aborted", file=sys.stderr)
        return

    if cwd.find("ultracam") > -1 or cwd.find("ultraspec") > -1:
        itype = 'U'
    elif cwd.find("hipercam") > -1:
        itype = 'H'
    else:
        print("** digest: failed to find one of hipercam, ultracam or ultraspec in path")
        print("digest aborted", file=sys.stderr)
        return

    # regular expression to match run directory, e.g. 2018-10 [2018, October],
    # 2020-21 [season 2020 to 2021], 2021-P106 [ESO period 106], 2021A [GTC semester],
    # 'Others' [TNT other people]
    rdre = re.compile("^(Others|\d\d\d\d-(\d\d|P\d\d\d)|\d\d\d\d[AB])$")

    # regular expression to match input night directory, e.g. 2018_10_30
    ndre = re.compile("^\d\d\d\d_\d\d_\d\d$")

    # regular expression to match hipercam raw data file
    rhre = re.compile("^run\d\d\d\d\.fits$")

    # regular expression to match ultra(cam|spec) raw data file
    rure = re.compile("^run\d\d\d\.dat$")

    # get the run directories
    rdirs = [
        rdir
        for rdir in os.listdir('.')
        if rdre.match(rdir) and os.path.isdir(rdir)
    ]
    rdirs.sort()

    if not len(rdirs):
        print("\n** Found no run directories", file=sys.stderr)
        print("digest aborted", file=sys.stderr)
        return

    print("\nChecking run-by-run ...\n")

    for rdir in rdirs:

        print(f"Run directory = {rdir}")
        telescope = os.path.join(rdir, 'telescope')
        if not os.path.exists(telescope):
            print(f'ERROR: file called "telescope" not found in run directory = {run}')
            return

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

        # Try to run all basic checks before actually doing anything
        for ndir in ndirs:

            # run elementary checks that are fast
            night = os.path.basename(ndir)

            nnew = os.path.basename(ndir).replace('_','-')
            if os.path.exists(nnew):
                print(f'A directory called {nnew} already exists clashing with {ndir}')
                print('digest aborted',file=sys.stderr)
                return

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


            telescope = os.path.join(ndir, 'telescope')
            if os.path.exists(telescope):
                print(f'ERROR: file called "telescope" already exists in {ndir}')
                return

            # compile list of runs in the directory
            if itype == 'H':
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


            if not skip_md5:

                md5 = os.path.join(ndir, f"MD5SUM_{night}")
                if os.path.isfile(md5):
                    print(f"... found the md5sum file = {md5}")

                    # extract runs from the MD5SUM file
                    mruns = []
                    with open(md5) as fin:
                        for line in fin:
                            hash, name = line.split()
                            if name.endswith(".xml") or name.endswith(".dat") or name.endswith(".fits"):
                                mruns.append(name[: name.rfind(".")])

                else:
                    mruns = None
                    print(
                        f"** no md5sum file called {md5} found", file=sys.stderr
                    )
                    if strict:
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
                    return


        if not skip_md5:

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

        print('\nChecks passed; will now carry out the work.')

        for ndir in ndirs:
            night = os.path.basename(ndir)
            nnight = night.replace("_", "-")
            rdir = os.path.dirname(ndir)

            if os.path.exists(nnight):
                print(f'{nnight} already exists and will not be over-written')
                continue

            # move the night directory up a level _ to -
            os.rename(ndir, nnight)
            print(f"\nmv {ndir} {nnight}")

            # link back into the run directory
            ltarget = os.path.join('..',nnight)
            lname = os.path.join(rdir, nnight)
            os.symlink(ltarget, lname)
            print(f"ln -s {ltarget} {lname}")

            # make a link in each night to the telescope
            # file in the run directory
            ltarget = os.path.join('..',rdir,'telescope')
            lname = os.path.join(nnight, 'telescope')
            os.symlink(ltarget, lname)
            print(f"ln -s {ltarget} {lname}")

            # copy the hand-written log file to a new version without
            # underscores or 'log' in the name
            oldlog = os.path.join(nnight, f"{night}_log.dat")
            newlog = os.path.join(nnight, f"{nnight}.dat")
            if os.path.exists(oldlog):
                if os.path.exists(newlog):
                    print(f"New style hand log {newlog} exists and will not be over-written")
                else:
                    shutil.copyfile(oldlog, newlog)
                    print(f"cp {oldlog} {newlog}")
            else:
                print(f"No hand log {oldlog} found")

    print('\nAll done!')
