"""
Exception classes
"""

__all__ = ['UltracamError','UendError','PowerOnOffError']

class UltracamError(Exception):
    """For throwing exceptions from the ultracam module"""
    pass

class UendError(UltracamError):
    """
    Exception for the standard way to reach the end of a data
    file (failure to read the timing bytes). This allows the
    iterator to die silently in this case while  raising
    exceptions for less anticipated cases.
    """
    pass

class PowerOnOffError(UltracamError):
    """
    Exception for trying to read a power on/off
    """
    pass
