"""
Classes and functions of general use
"""

import math
from .core import *

__all__ = ('Vec2D', 'add_extension', 'sub_extension', 'rgb')

class Vec2D:
    """Simple 2D vector class."""

    def __init__(self, x=0., y=0.):
        self.x = x
        self.y = y

    def length(self):
        """Returns the Euclidean length"""
        return math.sqrt(self.x**2+self.y**2)

    def unit(self):
        """Returns a unit vector version of the vector"""
        dist = self.length()
        if dist > 0:
            x = self.x / dist
            y = self.y / dist
            return Vec2D(x,y)
        else:
            return ValueError('cannot normalise a zero vector')

    def __iadd__(self, vec):
        """+=: in place addition"""
        self.x += vec.x
        self.y += vec.y
        return self

    def __add__(self, vec):
        """+: addition"""
        return Vec2D(self.x+vec.x, self.y+vec.y)

    def __isub__(self, vec):
        """-=: in place subtraction"""
        self.x -= vec.x
        self.y -= vec.y
        return self

    def __sub__(self, vec):
        """-: subtraction"""
        return Vec2D(self.x-vec.x, self.y-vec.y)

    def __imul__(self, const):
        """*=: in place multiplication by a constant"""
        self.x *= const
        self.y *= const
        return self

    def __mul__(self, const):
        """*: post multiplication by a constant"""
        return Vec2D(const*self.x, const*self.y)

    def __rmul__(self, const):
        """*: pre multiplication by a constant"""
        return Vec2D(const*self.x, const*self.y)

    def __itruediv__(self, const):
        """Divides by a constant"""
        self.x /= const
        self.y /= const

    def __repr__(self):
        return 'Vec2D(x={:f}, y={:f})'.format(self.x, self.y)

def add_extension(fname, ext):
    """Add extension ext to a file name if it is not already there, and returns
    the revised name

    """
    if len(ext) and not fname.endswith(ext):
        return '{}{}'.format(fname, ext)
    else:
        return fname

def sub_extension(fname, ext):
    """Subtracts extension ext from a file name if it is present, and returns
    the revised name

    """
    if fname.endswith(ext):
        return fname[:-len(ext)]
    else:
        return fname

def rgb(cname):
    """Returns the RGB tuple associated with colour name 'cname', following
    the colour definitions used by 'reduce'.

    Returns a ValueError if cname is not recognised.
    """
    if cname not in CNAMS:
        raise ValueError('colour = {:s} not recognised'.format(cname))
    else:
        return CIS[CNAMS[cnam]]


def script_args(args):
    """
    This is a small helper method that is used at the start of most
    of the HiPERCAM entry point scripts such as 'reduce' and 'grab'.

    Generally these are called with one argument 'args' set to None, in which
    case this takes the arguments from the system (which start with the
    command name), extract the root of the command name, and then return with
    that and the rest of the arguments. If args is not None, then it is
    expected that the first argument is the command name. 
    """
    if args is None:
        args = sys.argv
    command = os.path.split(args.pop(0))[1]
    return (command, args)
