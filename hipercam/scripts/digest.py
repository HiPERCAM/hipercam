import sys
import os
import shutil
import time
import glob
import re
import subprocess
import getpass 

import numpy as np
from astropy.time import Time, TimeDelta

import hipercam as hcam
from hipercam import cline, utils, spooler
from hipercam.cline import Cline

__all__ = ['digest',]

#######################################################
#
# digest -- ingests hipercam data for archival purposes
#
#######################################################

def digest(args=None):
    """Ingests hipercam data from the telescope for archiving purposes.

    ThiIt should be of little interest to most users. It does what the scripts
    checker, import_data, make_log_dirs and make_derived_dirs do for ULTRACAM
    data, but all in one go to reduce the need for thought.

    HiPERCAM run data from the telescope comes in a standarized form. The
    purpose of this script is run a few checks and set up a different
    standardised directory structure for it. It must be run inside a directory
    'hipercam' which should contain a sub-directory 'raw_data' which contains
    run directories of the form 2017-10 which themselves contain
    night-by-night directories of the form YYYY_MM_DD. In these, it expects to
    find files in these of the form 'run1234.fits', (always 4 digits), a file
    called MD5SUM_YYYY_MM_DD (i.e. matching the directory date), and a file
    called YYYY_MM_DD_log.dat. It will stop if any of these are not found and
    it is then up to the user to do something about it. It will also check
    that every run file appears in the MD5SUM file. Assuming these initial
    checks are past, it will then run 'md5sum -c' to check that the file
    md5sums match.

    The script has no arguments; its purpose is do a series of mundane and
    irritating tasks with minimal input from the user. It stops on problems
    and tries to report in enough detail to help them to get fixed.
    """

username = getpass. getuser() 
    RAW = 'raw_data'
    DERIVED = 'derived_data'
    LOGS = 'logs'

    # check that we are in the right place
    subdirs = [sdir for sdir in os.listdir() if os.path.isdir(sdir)]

    if os.path.basename(os.getcwd()) != 'hipercam' or RAW not in subdirs \
            or DERIVED not in subdirs or LOGS not in subdirs:
        print("""
** digest must be run in a directory called 'hipercam' which has sub-directories
'{:s}', '{:s}' and '{:s}'
""".format(RAW,DERIVED,LOGS),file=sys.stderr)
        print('digest aborted',file=sys.stderr)
        return

    # regular expressions to match runs, nights & files
    rdre = re.compile('^\d\d\d\d-\d\d$')
    ndre = re.compile('^\d\d\d\d_\d\d_\d\d$')
    rre = re.compile('^run\d\d\d\d\.fits$')

    # get the run directories
    rdirs = [os.path.join(RAW, rdir) for rdir in os .listdir(RAW) if \
                 os.path.isdir(os.path.join(RAW, rdir))]
    rdirs.sort()

    if len(rdirs):
        print('\nFound the following run directories: {!s}'.format(', '.join(rdirs)))
    else:
        print('\n** Found no run directories',file=sys.stderr)
        print('digest aborted',file=sys.stderr)
        return

    print('\nChecking run-by-run ...\n')

    for rdir in rdirs:
        print('Run directory = {:s}'.format(rdir))

        # loop through the run directories

        # look for a file called 'telescope'. They should all have one.
        telescope = os.path.join(rdir, 'telescope')
        if os.path.isfile(telescope):
            print('Found file = {:s}'.format(telescope))
        else:
            print('** Could not find file = {:s}'.format(telescope),file=sys.stderr)
            print('This should have the name of the telescope inside it',file=sys.stderr)
            print('digest aborted',file=sys.stderr)
            return


        # look for the night-by-night directories. If not present in matching form,
        # skip to the next on the basis that it has been done already
        ndirs = [os.path.join(rdir, ndir) for ndir in os.listdir(rdir) if \
                     ndre.match(ndir) and \
                     os.path.isdir(os.path.join(rdir, ndir))]
        ndirs.sort()

        if len(ndirs):
            tdirs = [os.path.basename(ndir) for ndir in ndirs]
            print('\nFound night directories {:s} in {:s}'.format(', '.join(tdirs),rdir))
        else:
            print('\nFound no night directories in '
                  '{:s}; will assume it has already been processed'.format(rdir))


        # match runs of specific form 'run1234.fits' (4 digits)
        print('\nChecking night-by-night ...\n')

        # initial checks
        for ndir in ndirs:

            print('Night directory = {:s}'.format(ndir))
            night = os.path.basename(ndir)
            log = os.path.join(ndir, '{:s}_log.dat'.format(night))
            if os.path.isfile(log):
                print('... found the log file = {:s}'.format(log))
            else:
                print('** no log file called {:s} found'.format(log), file=sys.stderr)
                print('digest aborted',file=sys.stderr)
                return

            md5 = os.path.join(ndir, 'MD5SUM_{:s}'.format(night))
            if os.path.isfile(md5):
                print('... found the md5sum file = {:s}'.format(md5))
            else:
                print('** no md5sum file called {:s} found'.format(md5), file=sys.stderr)
                print('digest aborted',file=sys.stderr)
                return

            # compile list of runs in the directory
            runs = [run[:-5] for run in os.listdir(ndir) if \
                        os.path.isfile(os.path.join(ndir, run)) and \
                        rre.match(run)]

            # extract runs from the log file
            lruns = []
            with open(log) as fin:
                for line in fin:
                    name = line.split()[0]
                    lruns.append(name)

            if set(runs) != set(lruns):
                print('** the runs in the log file do not match those in the directory',file=sys.stderr)
                print('Runs not in common are: {!s}'.format(set(runs) ^ set(lruns)))
                print('digest aborted',file=sys.stderr)
                return

            print('... found all the runs listed in the log file')

            # extract runs from the MD5SUM file
            mruns = []
            with open(md5) as fin:
                for line in fin:
                    hash, name = line.split()
                    mruns.append(name[:-5])

            if set(runs) != set(mruns):
                print('The runs in the md5sum file do not match those in the directory',file=sys.stderr)
                print('Runs not in common are: {!s}'.format(set(runs) ^ set(mruns)))
                print('digest aborted',file=sys.stderr)
                return

            print('... found all the runs listed in the md5sum file')

    # Now the md5sum checks
    print('\nNow running md5sum checks; these might take a while.\n')

    for rdir in rdirs:
        print('Run directory = {:s}'.format(rdir))

        # look for the night-by-night directories. If not present in matching form,
        # skip to the next on the basis that it has been done already
        ndirs = [os.path.join(rdir, ndir) for ndir in os.listdir(rdir) if \
                     ndre.match(ndir) and \
                     os.path.isdir(os.path.join(rdir, ndir))]
        ndirs.sort()

        for ndir in ndirs:
            print('Night directory = {:s}'.format(ndir))
            night = os.path.basename(ndir)

            """
            md5 = os.path.join(ndir, 'MD5SUM_{:s}'.format(night))
            with open(md5) as fin:
                for line in fin:
                    hash, run = line.split()
                    fname = os.path.join(ndir, run)

                    # run the md5sum
                    output = subprocess.check_output(['md5sum', fname]).decode()
                    hashc = output.split()[0]
                    if hashc == hash:
                        print('md5sum of {:s} is OK'.format(fname))
                    else:
                        print('** md5sum of {:s} is NOT OK!'.format(fname),file=sys.stderr)
                        print('digest aborted')
                        return
            """

    # Now change the data structure into my standard form:
    #
    # 1) Rename the night directories, replacing underscores in directory with '-' signs.
    # 2) Move each one up one level
    # 2) Create a soft-link to it in the night directory
    # 4) Copy the hand log of form yyyy_mm_dd_log.dat to yyyy-mm-dd.dat

    for rdir in rdirs:
        print('\n\nRun directory = {:s}'.format(rdir))

        # look for the night-by-night directories. If not present in matching form,
        # skip to the next on the basis that it has been done already
        ndirs = [os.path.join(rdir, ndir) for ndir in os.listdir(rdir) if \
                     ndre.match(ndir) and \
                     os.path.isdir(os.path.join(rdir, ndir))]
        ndirs.sort()

        for ndir in ndirs:
            night = os.path.basename(ndir)
            nnight = night.replace('_','-')
            nndir = os.path.join(RAW, nnight)
            lndir = os.path.join(os.path.dirname(ndir),
                                 nnight)

            # move the data directory up a level _ to -
#            os.rename(ndir, nndir)
            print('\nmv {:s} {:s}'.format(ndir, nndir))

            # link back into the run directory
#            os.symlink(nndir,lndir)
            print('ln -s {:s} {:s}'.format(nndir, lndir))

            # copy the hand-written log file to a new version without
            # underscores or 'log' in the name
            oldlog = os.path.join(nndir, '{:s}_log.dat'.format(night))
            newlog = os.path.join(nndir, '{:s}.dat'.format(nnight))
#            shutil.copyfile(oldlog, newlog)
            print('cp {:s} {:s}'.format(oldlog, newlog))

            # make a directory for log stuff
            logndir = os.path.join(LOGS, nnight)
#            os.mkdir(logndir)
            print('mkdir {:s}'.format(logndir))

            # create a link to the hand-written log in the log directory
            lname = '{:s}.dat'.format(nnight)
            file = os.path.join('..','..',nndir,lname)
            link = os.path.join(logndir,lname)
#            os.symlink(file,link)
            print('ln -s {:s} {:s}'.format(file, link))

            # create a link to the data directory in the log directory
            file = os.path.join('..','..',nndir)
            link = os.path.join(logndir,'data')
#            os.symlink(file,link)
            print('ln -s {:s} {:s}'.format(file, link))


            # make a directory for derived data
            derndir = os.path.join(DERIVED, nnight)
#            os.mkdir(derndir)
            print('mkdir {:s}'.format(derndir))

            # create a link to the data directory in the derived directory
            file = os.path.join('..','..',nndir)
            link = os.path.join(derndir,'data')
#            os.symlink(file,link)
            print('ln -s {:s} {:s}'.format(file, link))


            
