#!/usr/bin/env python

"""
Filters HiPERCAM logs for runs suited to making fringe maps. Basically requires full frames,
longish runs that were (probably) nodded, no binning. This works on the spreadsheet hipercam
log produced by 'hlogger' called 'hipercam-log.xlsx'. One line of selection criteria should
be easily adaptable.
"""

import numpy as np
import pandas as pd
from hipercam.utils import format_hlogger_table

if __name__ == '__main__':

    log = pd.read_excel(
        'hipercam-log.xlsx',
        dtype=object
    )

    # Filter the logs. Criteria:
    #
    # -- at least 5 frames
    # -- Full frame readout
    # -- Nodding mode enabled
    # -- Slide position suitable
    # -- run type = data
    # -- longish exposures
    # -- 1x1 binning

    flog = log.loc[
        (log['Nframe'] >= 5) &
        (log['Readout\nmode'] == 'FULL') &
        (log['Nod'] == 'On') &
        ((log['FPslide'] > 1000) | (log['FPslide'] == -99)) &
        (log['Run\ntype'] == 'data') &
        (log['Cadence\n(sec)'] > 10) &
        (log['XxY\nbin'] == '1x1') &
        (log['Sun alt\nstart'] < -18) &
        (log['Sun alt\nend'] < -18) &
        ((log['Moon\nphase'] < 0.3) | ((log['Moon alt\nstart'] < 0) & (log['Moon alt\nend'] < 0)))
    ]

    format_hlogger_table("hipercam-log-fringe-runs.xlsx", flog)

