import sys
import os
from time import gmtime, strftime

import hipercam as hcam
from hipercam import cline, utils
from hipercam.cline import Cline

# get hipercam version to write into the reduce file
from pkg_resources import get_distribution, DistributionNotFound
try:
    hipercam_version = get_distribution('hipercam').version
except DistributionNotFound:
    hipercam_version = 'not found'

__all__ = ['rupdate',]

################################################
#
# rupdate -- updates an old reduce file
#
################################################

def rupdate(args=None):
    """``rupdate rfile``

    As changes are made to reduce, the reduce files can become obsolete. This
    script tries to bring old reduce files up to date by adding in the new
    options but in a way that should give the old behaviour. The file is modified
    in place.

    Parameters:

        rfile    : string
           the output reduce file created using |setaper|. This will be read
           for the targets. The main target will be assumed to have been
           called '1', the main comparison '2'. If there is a '3' it will be
           plotted relative to '2'; all others will be ignored for plotting
           purposes.

    """
    command, args = utils.script_args(args)

    with Cline('HIPERCAM_ENV', '.hipercam', command, args) as cl:

        # register parameters
        cl.register('rfile', Cline.GLOBAL, Cline.PROMPT)

        # get inputs

        # the reduce file
        rfile = cl.get_value(
            'rfile', 'reduce output file',
            cline.Fname('reduce.red',hcam.RED,cline.Fname.OLD)
        )

    lines = []
    fversion = False
    nversion = 0
    while open(rfile) as fin:
        for line in fin:
            if line.startswith('version ='):
                version = line[9:].strip()
                if version == hcam.REDUCE_FILE_VERSION:
                    print('reduce file = {:s} is up to date'.format(file))
                    exit(0)
                elif version == '20181107':
                    # update version number
                    line = 'version = {:s}\n'.format(hcam.REDUCE_FILE_VERSION)

                    # record version
                    nversion = 1

            lines.append(line)

    # Write out modified file
    while open(rfile,'w') as fout:
        for line in lines:
            fout.write(line)

        fout.write("""
# The next section was automatically added by rupdate to update
# from an old version of reduce file.
""")

        if nversion == 1:
            fout.write('skipbadt = no\n")

    print('Updated reduce file = {:s}'.format(rfile))
