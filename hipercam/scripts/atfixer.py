import sys
import os
import re
import time
import hipercam as hcam

__all__ = ['atfixer',]

#########################################################
#
# atfixer -- runs tfixer on a series of night directories
#
#########################################################

def atfixer(args=None):
    """``atfixer``

    Specialist script to look for timing problems in all runs in 
    night directories.

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

        # change working directory to the tbytes one
        os.chdir(ndir)

        # read all the runs
        fnames = [fname for fname in os.listdir('.') if rnre.match(fname)]
        fnames.sort()

        for fname in fnames:
            if fname.endswith('.fits'):
                source = 'hl'
                run = fname [:-5]
            else:
                source = 'ul'
                run = fname [:-4]

            args = [None, 'prompt', source, run, '6', 'no', 'yes']

            try:
                hcam.scripts.tfixer(args)
            except hcam.HipercamError as err:
                print('ERROR =',err)

        print('Finished',ndir)

        # change up the working directory
        os.chdir('..')
