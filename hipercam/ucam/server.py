"""
ATC server access code extra. See Raw.py for more.
"""

import os
import urllib.request

__all__ = ['get_nframe_from_server', 'URL',]

# The ATC FileServer recognises various GET requests (look for 'action=' in
# the code) which are accessed using urllib. The following code is to allow
# urllib to connect to a local version of the FileServer which does not work
# without this code.

# prevent auto-detection of proxy settings
proxy_support = urllib.request.ProxyHandler({})
opener = urllib.request.build_opener(proxy_support)
urllib.request.install_opener(opener)

# Set the URL of the FileServer from the environment or default to localhost.
URL = os.environ['ULTRACAM_DEFAULT_URL'] if 'ULTRACAM_DEFAULT_URL' in os.environ else \
      'http://localhost/'

def get_nframe_from_server(run):
    """
    Returns the number of frames in the run via the FileServer

    Argument::

       run : (string)
          ULTRACAM run name as in 'run003'
    """

    # Get from FileServer
    full_url = URL + run + '?action=get_num_frames'
    resp = urllib.request.urlopen(full_url).read()

    # Parse the response
    loc = resp.find('nframes="')
    if loc > -1:
        end = resp[loc+9:].find('"')
        return int(resp[loc+9:loc+9+end])
    else:
        raise ValueError('get_nframe_from_server: failed to parse server response to ' + full_url)

def get_runs_from_server(dir=None):
    """
    Returns with a list of runs from the server

    Argument::

       dir : (string)
          name of sub-directory on server
    """

    # get from FileServer
    if dir is None:
        full_url = URL + '?action=dir'
    else:
        full_url = URL + dir + '?action=dir'
    resp = urllib.request.urlopen(full_url).read()

    # parse response from server
    ldir = resp.split('<li>')
    runs = [entry[entry.find('>run')+1:entry.find('>run')+7] for entry in ldir
            if entry.find('getdata">run') > -1]
    runs.sort()
    return runs
