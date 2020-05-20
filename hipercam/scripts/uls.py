import os
import re
import requests

__all__ = [
    "uls",
]

###################################################################
#
# uls -- 'ls'-like listing of runs available on the ULTRACAM server
#
###################################################################


def uls(args=None):
    """Gives an 'ls'-like listing of the runs available on the ULTRACAM file
    server. Just invoke as 'uls' with no arguments.

    See |hls| for the |hiper| equivalent.
    """

    # if ever want to get sub-directory list:
    #
    #    if dir is None:
    #        full_url = URL + '?action=dir'
    #    else:
    #        full_url = URL + dir + '?action=dir'

    url = (
        os.environ.get("ULTRACAM_DEFAULT_URL", "http://localhost:8007/") + "?action=dir"
    )

    # Regular expression to pick out file name from the server response
    fre = re.compile(r">(run\d\d\d)<")

    # send off the url
    response = requests.get(url, timeout=1)

    # splits up response, searches for file names
    fnames = []
    for line in response.text.split("\n"):
        m = fre.search(line)
        if m:
            fnames.append(m.group(1))

    # sort
    fnames.sort()

    # output to terminal
    if len(fnames):
        print("\n".join(fnames))
    else:
        print('No runs returned by the ULTRACAM server; perhaps you need to use "hls"?')
