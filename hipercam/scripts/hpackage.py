import sys
import os
import signal
import re
import shutil
import subprocess

import numpy as np
from astropy.time import Time
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

import hipercam as hcam
from hipercam import cline, utils, spooler, defect, fringe
from hipercam.cline import Cline

__all__ = [
    "hpackage",
]

##############################################################
#
# hpackage -- bundles up standard data from a run on an object
#
##############################################################


def hpackage(args=None):
    """``hpackage runs label``

    ``hpackage`` looks for standard reduce data products and bundles
    them up for users. The idea is to copy all the files needed to be
    bale to re-run the reduction with the pipeline, while also adding
    a few helpfule extras. Given 'run123' for instance, it looks for:

      run123.hcm -- typically the result from a run of |averun|
      run123.ape -- file of photometric apertures
      run123.red -- reduce file as made by |genred|
      run123.log -- result from |reduce|

    It also looks for calibration files inside the reduce file and
    copies them. It requires them to be within the same direcory and
    will fail if they are not.

    It produces several extra files which are:

      run123.fits -- FITS version of the log file
      run123_ccd1.fits -- joined-up ds9-able version of run123.hcm
                          (and ccd2 etc)
      run123_ccd1.reg -- ds9-region file representing the apertures
                         from run123.ape

    Arguments:

      run : str
         Series of run names of the ones to copy, separated by spaces.

      label : str
         Name for the directory to store all files forming root of
         output tar file

    """

    command, args = utils.script_args(args)
    FEXTS = (hcam.HCAM, hcam.APER, hcam.LOG, hcam.RED)

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("runs", Cline.LOCAL, Cline.PROMPT)
        cl.register("label", Cline.LOCAL, Cline.PROMPT)

        runs = cl.get_value(
            "runs", "run names [space separated]",
            'run005'
        )
        runs = runs.split()
        for run in runs:
            if os.path.dirname(run) != '':
                raise hcam.HipercamError(
                    'hpackage only runs on files in the working directory'
                )
            for fext in FEXTS:
                if not os.path.exists(run + fext):
                    raise hcam.HipercamError(
                        f'could not find {run+fext}'
                    )

        label = cl.get_value(
            "label", "label for output tar file",
            'hdata'
        )

    # Make name of temporary directory
    tdir = utils.temp_dir()
    tmpdir = os.path.join(tdir, label)

    with CleanUp(tmpdir) as cleanup:

        # create directory
        os.makedirs(tmpdir, exist_ok=True)
        print(f'Will write files to {tmpdir}')

        # Get on
        for run in runs:

            # strip extension
            root = os.path.splitext(run)[0]

            # need to read the file to determine
            # the number of CCDs
            print(
                run,root,utils.add_extension(run,hcam.HCAM)
            )
            mccd = hcam.MCCD.read(
                utils.add_extension(run,hcam.HCAM)
            )

            # convert the  hcm and ape files using joinup
            args = [
                None,'prompt','list','hf',run,'no'
            ]
            if len(mccd) > 1:
                args += ['0']
            args += [
                root + hcam.APER,
                'none','none','none','none','no',
                'float32',str(100),str(100),'no','rice',
                tmpdir
            ]
            hcam.scripts.joinup(args)

            # convert log to fits as well
            args = [
                None,'prompt',run,'h',tmpdir
            ]
            hcam.scripts.hlog2fits(args)

            # copy standard files over
            for fext in FEXTS:
                source = utils.add_extension(root,fext)
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

            # now the calibrations
            rfile = hcam.reduction.Rfile.read(run + hcam.RED)
            if rfile.bias is not None:
                source = utils.add_extension(
                    rfile.bias.head['FILENAME'], hcam.HCAM
                )
                if os.path.dirname(source) != '':
                    raise HipercamError(
                        f'bias = {source} is not in the present working directory'
                    )
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

            if rfile.dark is not None:
                source = utils.add_extension(
                    rfile.dark.head['FILENAME'], hcam.HCAM
                )
                if os.path.dirname(source) != '':
                    raise HipercamError(
                        f'dark = {source} is not in the present working directory'
                    )
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

            if rfile.flat is not None:
                source = utils.add_extension(
                    rfile.flat.head['FILENAME'], hcam.HCAM
                )
                if os.path.dirname(source) != '':
                    raise HipercamError(
                        f'flat = {source} is not in the present working directory'
                    )
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

            if rfile.fmap is not None:

                source = utils.add_extension(
                    rfile.fmap.head['FILENAME'], hcam.HCAM
                )
                if os.path.dirname(source) != '':
                    raise HipercamError(
                        f'fringe map = {source} is not in the present working directory'
                    )
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

                if rfile.fpair is not None:

                    source = utils.add_extension(
                        rfile.fpair.head['FILENAME'], hcam.HCAM
                    )
                    if os.path.dirname(source) != '':
                        raise HipercamError(
                            f'fringe peak/trough pair file = {source} is not in the present working directory'
                        )
                    target = os.path.join(tmpdir,source)
                    shutil.copyfile(source, target)
                    print(f'copied {source} to {target}')

        # tar up the results
        args = ['tar','cvfz',f'{label}.tar.gz','-C',tdir,label]
        subprocess.run(args)

class CleanUp:
    """
    Context manager to handle temporary files
    """
    def __init__(self, tmpdir):
        self.tmpdir = tmpdir

    def _sigint_handler(self, signal_received, frame):
        print("\nhpackage aborted")
        sys.exit(1)

    def __enter__(self):
        signal.signal(signal.SIGINT, self._sigint_handler)

    def __exit__(self, type, value, traceback):
        print(f'removing temporary directory {self.tmpdir}')
        shutil.rmtree(self.tmpdir)
