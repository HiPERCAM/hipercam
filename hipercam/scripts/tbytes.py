import sys
import os

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

__all__ = ['tbytes',]

############################################################
#
# tbytes -- strips timing bytes out of a run, writes to disk
#
############################################################

def tbytes(args=None):
    """``tbytes [source] run``

    Reads all timing bytes from a |hiper|, ULTRACAM or ULTRASPEC run
    and dumps them to a disk file. Designed as a safety fallback for
    correcting timing issues where one wants to manipulate the raw data.

    Parameters:

        source : string [hidden]
           Data source, two options:

              | 'hl' : local HiPERCAM FITS file
              | 'ul' : ULTRA(CAM|SPEC) server

        run : string
           run number to access, e.g. 'run0034'. This will also
           be used to generate the name for the timing bytes file
           (extension '.tbts'). If a file of this name already
           exists, the script will abort. The timing bytes file
           will be written to the present working directory, not
           necessarily the location of the data file.

    .. Warning::

       This routine does not yet work with ULTRACAM data.

    """

    command, args = utils.script_args(args)

    # get the inputs
    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('source', Cline.GLOBAL, Cline.HIDE)
        cl.register('run', Cline.GLOBAL, Cline.PROMPT)

        # get inputs
        source = cl.get_value(
            'source', 'data source [hl, ul]',
            'hl', lvals=('hl','ul')
        )

        run = cl.get_value('run', 'run name', 'run005')

    ofile = os.path.basename(run) + hcam.TBTS
    if os.path.exists(ofile):
        print('File =',ofile,'already exists; tbytes aborted')
        return
    
    if source == 'hl':

        with hcam.hcam.Rtbytes(run) as rtbytes:
            with open(ofile,'wb') as fout:
                for tbytes in rtbytes:
                    fout.write(tbytes)

    elif source == 'ul':

        # Read the header .xml file
        rhead = hcam.ucam.Rhead(run)

        if rhead.isPonoff():
            print(run,'is a power on/off and will be ignored.')
            return

        ifile = run + '.dat'

        fsize = rhead.framesize
        ntbytes = 2*rhead.headerwords

        with open(ofile,'wb') as fout:
            with open(ifile,'rb') as fin:
                nframe = 1
                while True:
                    tbytes = fin.read(ntbytes)
                    if len(tbytes) != ntbytes:
                        print('Finished reading timing bytes')
                        break
                    fout.write(tbytes)

                    # step to start of next frame
                    fin.seek(fsize-ntbytes,1)
                    nframe += 1

        print('Found',nframe,'frames in',run,'\nWrote timing data to',ofile)
