import numpy as np
from trm import cline
from trm.cline import Cline

import hipercam as hcam

__all__ = ["checkflags"]


def checkflags(args=None):
    """``checkflags log``

    Summarise the flags stored in a |reduce| ASCII log for each CCD and
    aperture. Counts are per flag, so a measurement with multiple flags
    contributes to each count. Unrecognised bits are reported separately.
    The counts use the stored aperture flags, without adding flags inferred
    from other columns such as ``MJDok``.

    Parameters:

      log : str
         name of the |reduce| ASCII log file
    """
    command, args = cline.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:
        cl.register("log", Cline.LOCAL, Cline.PROMPT)
        log = cl.get_value(
            "log", "reduce log file to check", cline.Fname("hcam", hcam.LOG)
        )

    hlog = hcam.hlog.Hlog.rascii(log)
    flags = [(name, flag) for name, flag in hcam.FLAGS if flag in hcam.FLAG_MESSAGES]
    known_mask = 0
    for _, flag in flags:
        known_mask |= int(flag)

    for cnam, ccd in hlog.items():
        for apnam in hlog.apnames[cnam]:
            values = ccd[f"flag_{apnam}"]
            nflagged = np.count_nonzero(values)
            print(f"CCD {cnam}, aperture {apnam}: {nflagged}/{len(values)} flagged")
            for name, flag in flags:
                count = np.count_nonzero(values & flag)
                if count:
                    print(f"  {name}: {count} - {hcam.FLAG_MESSAGES[flag]}")

            present = int(np.bitwise_or.reduce(values, initial=0))
            unknown = present & (int(hcam.ANY_FLAG) ^ known_mask)
            while unknown:
                bit = unknown & -unknown
                count = np.count_nonzero(values & bit)
                print(f"  UNKNOWN 0x{bit:08x}: {count}")
                unknown ^= bit
