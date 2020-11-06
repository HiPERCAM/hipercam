#!/usr/bin/env python

"""
Filters HiPERCAM logs for runs suited to making fringe maps. Selection criteria:

  1) >= 5 frames in the run
  2) Full frame readout
  3) No binning
  4) Nodding / dithering enabled
  5) Slide out or undefined
  6) Run type == data
  7) Cadence > 10 secs
  8) No binning
  9) Sun < -18 at start and end
 10) Moon phase < 0.3 or the Moon is below the horizon.

This returns relatively few runs; manual work after that.
"""

import numpy as np
import pandas as pd
from hipercam.utils import format_hlogger_table

if __name__ == '__main__':

    log = pd.read_excel(
        'hipercam-log.xlsx',
        dtype=object
    )

    # Filter the logs.

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

