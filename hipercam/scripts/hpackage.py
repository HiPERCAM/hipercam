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

from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam import spooler, defect, fringe

__all__ = [
    "hpackage",
]

##############################################################
#
# hpackage -- bundles up standard data from a run on an object
#
##############################################################


def hpackage(args=None):
    """``hpackage runs dname tar``

    ``hpackage`` looks for standard reduce data products and bundles
    them up into a single directory and optionally creates a tar
    file. The idea is to copy all the files needed to be able to
    re-run the reduction with the pipeline, while also adding a few
    helpful extras where possible. This is to make export of a reduction
    relatively painless.

    Given 'run123' for instance, it looks for:

      run123.hcm -- typically the result from a run of |averun|. This *must*
                    exist for each run included.

      run123.ape -- file of photometric apertures. Also must exist.

      run123*.red -- reduce file as made by |genred| (* == any string)

      run123*.log -- result from |reduce| (* == any string)

    Note the '*' means that multiple .red / .log pairs are possible to
    allow for multiple alternative reductions. But they must start
    'run123' and end '.red' and '.log' to be included. There has to be
    at least one such pair. It also looks for calibration files inside
    the reduce file and copies them. It requires them to be within the
    same directory and will fail if they are not. The reason for this
    is that they get copied over to the same directory and the
    directory structure required to store them elsewhere is not
    created. I often use soft-links to point to calibration files for
    this reason.

    It produces several extra files which are:

      run123.fits -- FITS version of the log file

      run123_ccd1.fits -- joined-up ds9-able version of run123.hcm
                          (and ccd2 etc) [but only if the windows are
                          in sync.

      run123_ccd1.reg -- ds9-region file representing the apertures
                         from run123.ape [see above re synced windows]

      README -- a file of explanation.

    The files are all copied to a temporary directory.

    Arguments:

      run : str
         Series of run names of the ones to copy, separated by spaces.

      dname : str
         Name for the directory to store all files forming the root of
         any output tar file created from it.

      tar : bool
         Make a tar.gz file of the directory at the end; the directory and
         the files in it will be deleted. Otherwise, no tar file is made and
         the directory is left untouched. The directory will however be deleted
         if the program is aborted early.

    """

    command, args = cline.script_args(args)
    FEXTS = (hcam.HCAM, hcam.APER)

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("runs", Cline.LOCAL, Cline.PROMPT)
        cl.register("dname", Cline.LOCAL, Cline.PROMPT)
        cl.register("tar", Cline.LOCAL, Cline.PROMPT)

        runs = cl.get_value(
            "runs", "run names [space separated]",
            'run005'
        )
        runs = runs.split()

        # for storing of log / red pairs
        logreds = {}
        flist = os.listdir('.')

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

            # try to look for .red / .log pairs
            rered = re.compile(f'({run}.*)\{hcam.RED}$')
            relog = re.compile(f'({run}.*)\{hcam.LOG}$')
            reds, logs = [], []
            for fname in flist:
                m = rered.match(fname)
                if m:
                    reds.append(m.group(1))
                m = relog.match(fname)
                if m:
                    logs.append(m.group(1))

            logred = set(reds) & set(logs)
            if len(logred) == 0:
                raise hcam.HipercamError(
                    f'found no .log / .red pair starting "{run}"'
                )
            logreds[run] = logred

        dname = cl.get_value(
            "dname",
            "name of directory for storage of files (will be used for any tar file too)",
            'hdata'
        )

        tar = cl.get_value(
            "tar", "make a tar file (and delete temporary directory at end)?", True
        )

    # Make name of temporary directory
    tdir = utils.temp_dir()
    tmpdir = os.path.join(tdir, dname)

    with CleanUp(tmpdir, tar) as cleanup:

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
                run,root,cline.add_extension(run,hcam.HCAM)
            )
            mccd = hcam.MCCD.read(
                cline.add_extension(run,hcam.HCAM)
            )

            try:
                # convert the  hcm and ape files using joinup
                args = [
                    None,'prompt','list','hf',run,'no'
                ]
                if len(mccd) > 1:
                    args += ['0']
                args += [
                    root + hcam.APER,
                    'none','none','none','none','no',
                    'float32',str(100),str(100),
                    'no','rice',tmpdir
                ]
                hcam.scripts.joinup(args)
            except:
                print('failed to joinup hcm / aper files')

            # copy standard files over
            for fext in FEXTS:
                source = cline.add_extension(root,fext)
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

            # try to look for alternative .red / .log pairs
            for logred in logreds[run]:

                # copy .red and .log over
                source = cline.add_extension(logred,hcam.RED)
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

                source = cline.add_extension(logred,hcam.LOG)
                target = os.path.join(tmpdir,source)
                shutil.copyfile(source, target)
                print(f'copied {source} to {target}')

                # convert log to fits
                args = [
                    None,'prompt',logred,'h',tmpdir
                ]
                hcam.scripts.hlog2fits(args)

                # now the calibrations
                rfile = hcam.reduction.Rfile.read(logred + hcam.RED)
                csec = rfile['calibration']
                if rfile.bias is not None:
                    source = cline.add_extension(
                        csec['bias'], hcam.HCAM
                    )
                    if os.path.dirname(source) != '':
                        raise hcam.HipercamError(
                            f'bias = {source} is not in the present working directory'
                        )
                    target = os.path.join(tmpdir,source)
                    if os.path.exists(target):
                        print(f'copy skipped: {target} exists')
                    else:
                        shutil.copyfile(source, target)
                        print(f'copied {source} to {target}')

                if rfile.dark is not None:
                    source = cline.add_extension(
                        csec['dark'], hcam.HCAM
                    )
                    if os.path.dirname(source) != '':
                        raise hcam.HipercamError(
                            f'dark = {source} is not in the present working directory'
                        )
                    target = os.path.join(tmpdir,source)
                    if os.path.exists(target):
                        print(f'copy skipped: {target} exists')
                    else:
                        shutil.copyfile(source, target)
                        print(f'copied {source} to {target}')

                if rfile.flat is not None:
                    source = cline.add_extension(
                        csec['flat'], hcam.HCAM
                    )
                    if os.path.dirname(source) != '':
                        raise hcam.HipercamError(
                            f'flat = {source} is not in the present working directory'
                        )
                    target = os.path.join(tmpdir,source)
                    if os.path.exists(target):
                        print(f'copy skipped: {target} exists')
                    else:
                        shutil.copyfile(source, target)
                        print(f'copied {source} to {target}')

                if rfile.fmap is not None:
                    source = cline.add_extension(
                        csec['fmap'], hcam.HCAM
                    )
                    if os.path.dirname(source) != '':
                        raise hcam.HipercamError(
                            f'fringe map = {source} is not in the present working directory'
                        )
                    target = os.path.join(tmpdir,source)
                    if os.path.exists(target):
                        print(f'copy skipped: {target} exists')
                    else:
                        shutil.copyfile(source, target)
                        print(f'copied {source} to {target}')

                    if rfile.fpair is not None:
                        source = cline.add_extension(
                            csec['fpair'], hcam.FRNG
                        )
                        if os.path.dirname(source) != '':
                            raise hcam.HipercamError(
                                f'fringe peak/trough pair file = {source}'
                                ' is not in the present working directory'
                            )
                        target = os.path.join(tmpdir,source)
                        if os.path.exists(target):
                            print(f'copy skipped: {target} exists')
                        else:
                            shutil.copyfile(source, target)
                            print(f'copied {source} to {target}')

        readme = os.path.join(tmpdir,'README')
        with open(readme,'w') as fp:
            fp.write(README)

        # tar up the results
        args = ['tar','cvfz',f'{dname}.tar.gz','-C',tdir,dname]
        subprocess.run(args)

class CleanUp:
    """
    Context manager to handle temporary files
    """
    def __init__(self, tmpdir, tar):
        self.tmpdir = tmpdir
        self.delete = tar

    def _sigint_handler(self, signal_received, frame):
        print("\nhpackage aborted")
        self.delete = True
        sys.exit(1)

    def __enter__(self):
        signal.signal(signal.SIGINT, self._sigint_handler)

    def __exit__(self, type, value, traceback):
        if self.delete:
            print(f'removing temporary directory {self.tmpdir}')
            shutil.rmtree(self.tmpdir)
        else:
            print(f'All files have been written to {self.tmpdir}')

README = """
This tar file contains data products from the HiPERCAM pipeline, (but
might include HiPERCAM, ULTRACAM or ULTRASPEC data, depending on the
instrument in use). The aim is to provide all the files needed to
carry out a (re-)reduction with the pipeline command "reduce", along
with some support files that integrate with 'ds9' and 'fv'. For each
run, say "run123", you should find the following files:

  run123.ape -- JSON file of photometric apertures
  run123.hcm -- HiPERCAM file of CCD data (see note below)
  run123.log -- ASCII log from running "reduce", i.e. the photometry
  run123.red -- text file of reduction parameters
  run123.fits - a FITS-format version of run123.log which is usually
                easier to understand (e.g. look at it with 'fv').
  run123_ccd1.fits -- single joined up HDU version of CCD 1 from run123.hcm
  run123_ccd2.fits -- same for CCD 2, if there is one
  .
  .
  run123_ccd1.reg -- set of ds9 "regions" equivalent to run123.ape
  run123_ccd2.reg -- same for CCD 2, if there is one
  .
  .

although the _ccd1.fits & .reg files may not be present if the attempt
to joinup the data failed owing to non-synchronised windowed data.

You may also find some other "hcm" files with names like bias.hcm,
flat.hcm, dark.hcm, and fmap.hcm (which should come accompanied by
something like fpair.frng). These are calibration files that needed
by "reduce" (look inside "run123.red").

Notes:

1) "hcm" files are multi-HDU FITS files that can be looked at by ds9,
   but the joined up images should have a WCS in the case of HiPERCAM,
   and will be easier to start with.

2) Typically this package of files will come from the telescope where
   reduction is done on the fly. It is *extremely likely* that you can
   optimise the settings in the ".red" reduction file to improve your
   data. For instance the default sky annuli are usually quite large
   to accommodate possible seeing variations, but if in fact seeing is
   stable, they can probably be reduced in a re-reduction.

3) By convention, we assign aperture 1 to the main target, 2 to the
   brightest comparison. But note that the comparison may change between
   CCDs, so always look at all of them to be sure.

4) Given "run123_ccd2.fits" and "run123_ccd2.reg", the following command
   should display them with ds9, using the "zscale" intensity scaling
   option:

     ds9 run123_ccd2.fits -regions run123_ccd2.reg -zscale

5) The HiPERCAM pipeline can be obtained from:

       https://github.com/HiPERCAM/hipercam

If you encounter problems, please contact Tom Marsh.
"""
