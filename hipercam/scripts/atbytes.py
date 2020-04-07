import sys
import os
import re
import time
import hipercam as hcam

__all__ = ['atbytes',]

###########################################################################
#
# atbytes -- strips timing bytes out of all runs in a series of directories
#
###########################################################################

def atbytes(args=None):
    """``atbytes``

    Script to strip out and save all timing bytes of all runs in a
    series of directories of YYY-MM-DD form which should be
    sub-directories of the directory the script is run from. It will
    create (if necessary) a sub-directory of each of these call
    'tbytes' containing the timing bytes data. It can handle ULTRACAM,
    ULTRASPEC and HiPERCAM data.

    """

    # Specific formats for night directories and runs within them
    ndre = re.compile('^\d\d\d\d[_-]\d\d[_-]\d\d$')
    rnre = re.compile('^run\d\d\d\d\.fits|run\d\d\d\.xml$')

    # Find list of night directories
    ndirs = [dname for dname in os.listdir('.') \
             if ndre.match(dname) and os.path.isdir(dname)]
    ndirs.sort()

    for ndir in ndirs:

        print('\nStarted on',ndir)

        # create sub-directory where all timing info will be saved
        tbytes_dir = os.path.join(ndir,'tbytes')
        if not os.path.exists(tbytes_dir):
            os.mkdir(tbytes_dir)

        # read all the runs in the night directory
        fnames = [fname for fname in os.listdir(ndir) if rnre.match(fname)]
        fnames.sort()

        # change working directory to the tbytes one
        os.chdir(tbytes_dir)
        print('Changed working directory to',tbytes_dir)

        for fname in fnames:
            if fname.endswith('.fits'):
                source = 'hl'
                run = os.path.join('..',fname [:-5])
            else:
                source = 'ul'
                run = os.path.join('..',fname [:-4])
            tfile = run + hcam.TBTS
            if os.path.exists(tfile):
                print(tfile,'already exists; will not re-make')
            else:
                args = [None, 'prompt', source, run]
                try:
                    hcam.scripts.tbytes(args)
                except hcam.ucam.PowerOnOffError:
                    print('ignoring',run,'which is a Power On or Off')

        print('Finished',ndir)

        # change up the working directory to the tbytes one
        os.chdir('../..')
        print('Moved working directory up two levels')
