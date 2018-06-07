import os
import re
import requests

__all__ = ['hls',]

###################################################################
#
# hls -- 'ls'-like listing of runs available on the HiPERCAM server
#
###################################################################

def hls(args=None):
    """Gives an 'ls'-like listing of the runs available on the HiPERCAM file
    server. Just invoke as 'hls' with no arguments. This should work correctly
    whether it in fact its the HiPERCAM or ULTRACAM file server that is
    running.

    """

    url = 'http://' + \
          os.environ.get('HIPERCAM_DEFAULT_URL', 'http://localhost:8007/') + \
          '?action=dir'

    # Regular expression to pick out file name from the server response
    fre = re.compile(r'^run\d\d\d\d$')

    # send off the url
    response = requests.get(url, timeout=1)

    # splits up response, searches for file names
    fnames = [fname for fname in response.text.split('\n') if fre.search(fname)]

    # sort
    fnames.sort()

    # output to terminal
    if len(fnames):
        print('\n'.join(fnames))
    else:
        print('No runs returned by the HiPERCAM server; perhaps you need to use "uls"?')
