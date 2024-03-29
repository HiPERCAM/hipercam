import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from astropy.io import fits

from trm import cline
from trm.cline import Cline

import hipercam as hcam
from hipercam import core, spooler, defect

__all__ = [
    "ltrans",
]

###############################################
#
# ltrans -- computes transforms to align frames
#
###############################################


def ltrans(args=None):
    """``ltrans posns``

    Reads a file of target positions dumped by `ftargets` and uses it
    to derive the transformations needed to co-align images.

    Parameters:

        posns : string
           file of positions generated by `ftargets`

        cmax : float
           maximum number of counts at peak allowed for a target. Targets
           with greater than this will be ignored when computing the transform.

        emax : float
           maximum ratio a/b, major/minor axis elongation. Targets more elongated
           than this will be ignored.
    """

    command, args = cline.script_args(args)

    # get the inputs
    with Cline("HIPERCAM_ENV", ".hipercam", command, args) as cl:

        # register parameters
        cl.register("posns", Cline.GLOBAL, Cline.PROMPT)
        cl.register("cmax", Cline.GLOBAL, Cline.PROMPT)
        cl.register("emax", Cline.GLOBAL, Cline.PROMPT)

        # bias frame (if any)
        posns = cl.get_value(
            "posns",
            "file of target positions in a series of frames",
            cline.Fname("sources", hcam.SEP),
        )

        cmax = cl.get_value("cmax", "maximum peak counts for any target", 60000.0)

        emax = cl.get_value("emax", "maximum elongation (a/b)", 1.5, 1.0)

    with fits.open(posns) as hdul:
        for hdu in hdul[1:]:
            table = hdu.data

            nfs = table["nframe"]
            xs, ys = table["x"], table["y"]
            exs, eys = np.sqrt(table["errx2"]), np.sqrt(table["erry2"])

            fig = plt.figure()
            axes = fig.add_subplot(111)
            axes.set_aspect("equal")
            plt.errorbar(xs, ys, eys, exs, ".")
            plt.title('Raw positions from all frames')
            plt.show()

            nfmax = 0
            nfm = 0
            for nf in range(nfs.min(), nfs.max() + 1):
                ok = nf == nfs
                ntargs = len(nfs[ok])
                if ntargs > nfmax:
                    nfmax = ntargs
                    nfref = nf
                    okref = ok

            print(
                "Maximum number of targets = {:d} for frame {:d}".format(nfmax, nfref)
            )
            # reference X,Y values
            xrefs, yrefs = xs[okref], ys[okref]

            xts, yts = np.empty_like(xs), np.empty_like(ys)

            for nf in range(nfs.min(), nfs.max() + 1):
                ok = nf == nfs
                if nf == nfref:
                    xts[ok], yts[ok] = xs[ok], ys[ok]
                else:
                    theta, (xoff, yoff) = matchpos(xs[ok], ys[ok], xrefs, yrefs)
                    theta = np.radians(theta)
                    ct, st = np.cos(theta), np.sin(theta)
                    xts[ok] = ct * (xs[ok] - xoff) - st * (ys[ok] - yoff)
                    yts[ok] = st * (xs[ok] - xoff) + ct * (ys[ok] - yoff)

            fig = plt.figure()
            axes = fig.add_subplot(111)
            axes.set_aspect("equal")
            plt.errorbar(xts, yts, eys, exs, ".")
            plt.title('Corrected positions from all frames')
            plt.show()


def matchpos(x, y, xref, yref):
    """
    Match pairs of positions, looking for the nearest neighbour.

    Parameters
    ----------
    x, y: `numpy.ndarray`
        spots to match
    xref, yref: `numpy.ndarray`
        reference positions

    Returns
    -------
    theta, offset: `numpy.ndarray`
        rotation and offset that matches points to reference points
    """
    tree = cKDTree(np.column_stack((xref, yref)))
    seps, indices = tree.query(np.column_stack((x, y)))
    if np.all(seps < 0.01):
        # probably the same data
        mask = np.zeros_like(seps).astype("bool")
    else:
        sep_sigma = (seps - seps.mean()) / seps.std()
        mask = sep_sigma < 1.5

    return findBestRigidTransform(
        x[mask], y[mask], xref[indices][mask], yref[indices][mask]
    )


def findBestRigidTransform(x, y, xref, yref):
    """Given a set of points and *matching* reference points, find the
    optimal Rigid transform.

    A Rigid transform is a rotation and translation, i.e:

    B = R.A + t

    See https://scicomp.stackexchange.com/questions/6878/fitting-one-set-of-points-to-another-by-a-rigid-motion
    or http://nghiaho.com/?page_id=671

    for method.

    """

    A = np.column_stack((xref, yref))
    B = np.column_stack((x, y))

    # find centre of points
    centroidA = np.mean(A, axis=0)
    centroidB = np.mean(B, axis=0)

    # shift centres to origin
    Ashifted = A - centroidA
    Bshifted = B - centroidB

    # find dot product
    H = np.dot(Bshifted.T, Ashifted)

    # use SVD to find rotation
    u, s, v = np.linalg.svd(H)
    R = np.dot(v, u.T)
    if np.linalg.det(R) < 0:
        R = np.dot(v, np.array([1, -1]) * u.T)

    T = np.median(B - np.dot(A, R), axis=0)
    theta1 = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
    theta2 = np.degrees(-np.arctan2(R[0, 1], R[1, 1]))
    theta = 0.5 * (theta1 + theta2)
    return theta, T
