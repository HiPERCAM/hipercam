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
    """``atfixer mintim dcmax``

    Specialist script to look for timing problems in all runs in 
    night directories.

    Parameters:

        mintim : int
           Minimum number of frames to reuire a run to have before running tfixer on it

        dcmax : float
           Maximum differential in terms of cycle number from exact integers
           to use to indicate "failed" times. Pre 2010-02 needs to be looser
           than post 2010-02. Number of order 0.003 is appropriate.

    """

    if not sys.warnoptions:
        # to avoid some annoying warnings from utimer
        import warnings
        warnings.simplefilter("ignore")

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('mintim', Cline.LOCAL, Cline.PROMPT)
        cl.register('dcmax', Cline.LOCAL, Cline.PROMPT)

        mintim = cl.get_value('mintim', 'minimum number of times needed', 6, 4)
        dcmax = cl.get_value('dcmax', 'maximum cycle number offset from an integer', 0.003, 0)

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

            args = [None, 'prompt', source, run, str(mintim), str(dcmax), 'no', 'no', 'yes']

            try:
                hcam.scripts.tfixer(args)
            except hcam.HipercamError as err:
                print('ERROR =',err)

        print('Finished',ndir)

        # change up the working directory
        os.chdir('..')
