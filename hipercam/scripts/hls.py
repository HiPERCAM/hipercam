import os
import sys
import re
from urllib.request import urlopen

__all__ = ['hls',]

##########################################################
#
# hls -- 'ls'-like listing of runs available on the server
#
##########################################################

# URL of the FileServer
URL = 'http://{:s}'.format(os.environ['HIPERCAM_DEFAULT_URL']) \
      if 'HIPERCAM_DEFAULT_URL' in os.environ else \
         'http://localhost:8007/'

# need to pick up file names from either the hipercam or the ultracam
# server. The hipercam server simply returns a list of files in the form
# run####, one per line whereas the ultracam server surrounds it with html
# stuff.
fre = re.compile(r'^(run\d\d\d\d)$|>(run\d\d\d)<')

def hls(args=None):
    """Gives an 'ls'-like listing of the runs available on the HiPERCAM file
    server. Just invoke as 'hls' with no arguments. This should work correctly
    whether it in fact its the HiPERCAM or ULTRACAM file server that is
    running.

    """

    url = '{:s}?action=dir'.format(URL)

    with urlopen(url) as f:
        response = f.read().decode('utf-8')

    # splits up response, searches for file names of the
    # form runXXX[X], prints them as a list. If it finds
    fnames = []
    for line in response.split('\n'):
        m = fre.search(line)
        if m:
            fnames.append(m.group(0) if m.group(1) else m.group(2))

    # sort
    fnames.sort()

    # output to terminal
    if len(fnames):
        print('\n'.join(fnames))
    else:
        print('No runs returned by server')
