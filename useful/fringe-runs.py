#!/usr/bin/env python

"""
Filters HiPERCAM logs for runs suited to making fringe maps. Basically requires full frames,
longish runs that were (probably) nodded, no binning. This works on the spreadsheet hipercam
log produced by 'hlogger' called 'hipercam-log.xlsx'. One line of selection criteria should
be easily adaptable.
"""

import pandas as pd

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
        (log['XxY\nbin'] == '1x1')
    ]

    writer = pd.ExcelWriter("hipercam-log-fringe-runs.xlsx", engine='xlsxwriter')
    flog.to_excel(writer, sheet_name='Hipercam logs', index=False)

    # format
    worksheet = writer.sheets["Hipercam logs"]

    # Column widths. Should probably automate. Needs to track whatever is in 'hlogger'
    worksheet.set_column('A:A', 8)
    worksheet.set_column('B:B', 30)
    worksheet.set_column('C:D', 11)
    worksheet.set_column('F:G', 10)
    worksheet.set_column('H:I', 9)
    worksheet.set_column('J:J', 6)
    worksheet.set_column('M:N', 8)
    worksheet.set_column('O:P', 12)
    worksheet.set_column('Q:R', 8)
    worksheet.set_column('S:T', 18)
    worksheet.set_column('U:X', 4)
    worksheet.set_column('Y:Z', 5)
    worksheet.set_column('AA:AA', 4)
    worksheet.set_column('AB:AB', 7)
    worksheet.set_column('AC:AC', 5)
    worksheet.set_column('AD:AD', 6)
    worksheet.set_column('AE:AF', 7)
    worksheet.set_column('AG:AG', 26)
    worksheet.set_column('AH:AH', 15)
    worksheet.set_column('AI:AI', 14)
    worksheet.set_column('AJ:AJ', 18)

    # Fancy stuff with links in penultimate column
    cell_format = writer.book.add_format({'font_color': 'blue'})
    worksheet.set_column('AK:AK', 10, cell_format)

    # set links using a formula
    for nr in range(2,len(flog)+2):
        worksheet.write_formula(
            f'AK{nr}',
            f'=HYPERLINK("http://deneb.astro.warwick.ac.uk/phsaap/hipercam/logs/" & F{nr}, A{nr})'
        )

    worksheet.set_column('AL:AL', 200)

    worksheet.set_row(0,28)
    worksheet.freeze_panes(1, 0)
    worksheet.set_zoom(200)
    writer.save()


