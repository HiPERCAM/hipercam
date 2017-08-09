"""
Exception classes
"""

__all__ = ['HipercamError','HendError']

class HipercamError(Exception):
    """For throwing exceptions from the hipercam.hcam sub-module"""
    pass

class HendError(HipercamError):
    """
    Exception for the standard way to reach the end of a HiPERCAM raw data
    file (attempt to read out-of-range frame). This allows the iterator to die
    silently in this case while raising exceptions for less anticipated cases.
    """
    pass
