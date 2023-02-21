import sys
import os
from time import gmtime, strftime

from trm import cline
from trm.cline import Cline

import hipercam as hcam

# get hipercam version to write into the reduce file
from pkg_resources import get_distribution, DistributionNotFound

try:
    hipercam_version = get_distribution("hipercam").version
except DistributionNotFound:
    hipercam_version = "not found"

__all__ = [
    "rupdate",
]

################################################
#
# rupdate -- updates an old reduce file
#
################################################


def rupdate(args=None):
    """``rupdate rfile``

    As changes are made to 'reduce', old reduce files can become
    obsolete. This script tries to bring old reduce files up to date
    by adding in the new options but in a way that should give the old
    behaviour. The file is modified in place.

    There are some cases where a full update might not be possible and
    a further manual edit of the reduce file might be required.

    Parameters:

        rfile    : string
           the output reduce file created using |setaper|. This will be read
           for the targets. The main target will be assumed to have been
           called '1', the main comparison '2'. If there is a '3' it will be
           plotted relative to '2'; all others will be ignored for plotting
           purposes.

    """
    command, args = cline.script_args(args)

    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("rfile", Cline.GLOBAL, Cline.PROMPT)

        # get inputs

        # the reduce file
        rfile = cl.get_value(
            "rfile",
            "reduce output file",
            cline.Fname("reduce.red", hcam.RED, cline.Fname.OLD),
        )

    lines = []
    fversion = False
    nversion = 0
    with open(rfile) as fin:
        for line in fin:
            if line.startswith("version ="):
                # extract the version number
                version = line[9:].strip()
                if version.find(" ") > -1:
                    version = version[: version.find(" ")]
                if version.find("#") > -1:
                    version = version[: version.find("#")]

                if version == hcam.REDUCE_FILE_VERSION:
                    print(f"reduce file = {rfile} is up to date")
                    exit(0)

                elif version == "20181107":
                    # update version number
                    line = (
                        "version = {:s} # must be"
                        " compatible with the"
                        " version in reduce\n"
                    ).format(hcam.REDUCE_FILE_VERSION)

                    # Insert modified version and extra lines which go into the 'general'
                    # section along with the version number.
                    lines.append(line)
                    lines.append("\n# Next line was automatically added by rupdate\n")
                    lines.append("skipbadt = no\n\n")

                    # record version in case we need other actions later
                    nversion = 1

                elif version == "20200207":
                    # update version number
                    line = (
                        "version = {:s} # must be"
                        " compatible with the"
                        " version in reduce\n"
                    ).format(hcam.REDUCE_FILE_VERSION)
                    lines.append(line)

                    # record version in case we need other actions later
                    nversion = 2

                elif version == "20200223":
                    # update version number
                    line = (
                        "version = {:s} # must be"
                        " compatible with the"
                        " version in reduce\n"
                    ).format(hcam.REDUCE_FILE_VERSION)
                    lines.append(line)

                    # record version in case we need other actions later
                    nversion = 3

                elif version == "20200318":
                    # update version number
                    lines.append(
                        f"version = {hcam.REDUCE_FILE_VERSION} # must be"
                        " compatible with the version in reduce\n"
                    )

                    # record version in case we need other actions later
                    nversion = 4

                elif version == "20210523":
                    # update version number
                    lines.append(
                        f"version = {hcam.REDUCE_FILE_VERSION} # must be"
                        " compatible with the version in reduce\n"
                    )

                    # record version in case we need other actions later
                    nversion = 5

                elif version == "20210602":
                    # update version number
                    lines.append(
                        f"version = {hcam.REDUCE_FILE_VERSION} # must be"
                        " compatible with the version in reduce\n"
                    )

                    # record version in case we need other actions later
                    nversion = 6

                else:
                    print("Version = {:s} not recognised".format(version))
                    print("Aborting update; nothing changed.")
                    exit(1)

                if nversion <= 5:
                    lines.append("\n# Next lines were automatically added by rupdate\n")
                    lines.append("instrument = UNKNOWN # instrument-telescope\n")
                    lines.append("scale = UNKNOWN # scale, arcsec/unbinned pixel\n")

            elif line.startswith("fit_fwhm_min ="):
                lines.append(line)
                if nversion <= 3:
                    lines.append("\n# Next line was automatically added by rupdate\n")
                    lines.append(
                        "fit_fwhm_max = 1000 # Maximum FWHM, unbinned pixels\n"
                    )

            elif line.startswith("dark ="):
                lines.append(line)
                if nversion <= 4:
                    lines.append("\n# Next lines were automatically added by rupdate\n")
                    lines.append("fmap = # Fringe map, blank to ignore\n")
                    lines.append("fpair = # FringePair file, ignored if fringe blank\n")
                    lines.append("nhalf = 3 # Half-width, ignored if fringe blank\n")
                    lines.append("rmin = -2 # minimum ratio, ignored if fringe blank\n")
                    lines.append("rmax = 1 # maximum ratio ignored if fringe blank\n\n")

            elif line.startswith("scale =") and nversion <= 5:
                arr = line.split()
                try:
                    scale = float(arr[2])
                except:
                    scale = "UNKNOWN"

            elif line.startswith("[psf_photom]") and nversion <= 6:
                lines.append(line)
                lines.append("\n# Next 2 lines were automatically added by rupdate\n")
                lines.append("use_psf = no # use PSF photometry\n")
                lines.append(
                    "psf_model = 'moffat' # must be one of 'gaussian', 'moffat'\n"
                )

            else:
                # Default action is just to store save the line
                lines.append(line)

    # Write out modified file
    with open(rfile, "w") as fout:
        for line in lines:
            if line.startswith("scale =") and nversion <= 5 and scale != "UNKNOWN":
                fout.write(f"scale = {scale:.3f} # scale, arcsec per unbinned pixel\n")
            else:
                fout.write(line)

        # This could be the point at which extra lines are tacked
        # on to the end of the file.
        # if nversion == XX etc
        if nversion == 1 or nversion == 2:
            fout.write(
                """#
# Next lines were added automatically by 'rupdate' to bring this reduce file
# up to date. The changes are designed to produce the same behaviour as the
# old version of the reduce file as far as is possible.

[focal_mask]
demask = no
dthresh = 4
"""
            )

    print(f"Updated reduce file = {rfile}")
