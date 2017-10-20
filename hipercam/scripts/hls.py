import os
import sys
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

def hls(args=None):
    """Gives an 'ls'-like listing of the runs available on the HiPERCAM
    server. Just invoke as 'hls' with no arguments.
    """

    url = '{:s}?action=dir'.format(URL)

    with urlopen(url) as f:
        response = f.read().decode('utf-8')

    # splits up response, searches for file names of the
    # form runXXXX, prints them as a list.
    fnames = [fname for fname in response.split('\n') \
              if fname.startswith('run') and len(fname) == 7]

    # sort
    fnames.sort()

    # output to terminal
    print('\n'.join(fnames))



