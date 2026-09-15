"""Protocols and optional abstract base classes for supporting third-party
instruments as first-class citizens in the hipercam pipeline.

Plugin discovery is via the ``hipercam.instruments`` entry-point group.
Third-party packages register themselves with an entry pointing to a module
that satisfies :class:`InstrumentModuleProtocol`.

Example ``pyproject.toml`` registration::

    [project.entry-points."hipercam.instruments"]
    myinstrument = "mypackage.mymodule"

The registered module must expose:

* ``Rdata``     — data-reading class (iterable, callable, context-manager)
* ``supports_server_access`` — boolean indicating if server access is supported
* ``get_ccd_pars(resource, server=False)`` — returns CCD geometry as an
  :class:`~collections.OrderedDict`

Convenience base class  :class:`RdataBase` is provided for plugin authors 
who want to reduce boilerplate, but inheriting from it is **not** required.
"""
import warnings
import sys
from abc import ABC, abstractmethod
from typing import Optional

from hipercam import HipercamError

__all__ = [
    "IendError",
    "InstrumentModuleProtocol",
    "RdataBase",
    "discover_instruments",
    "available_instruments",
    "load_instrument",
]


# Shared end-of-data exception
class IendError(HipercamError):
    """Exception for the standard way to reach the end of a raw data file.

    Raised when an attempt is made to read a frame beyond the end of the
    available data. This allows iterators to die silently on normal
    end-of-file while still propagating unexpected errors.
    """

    pass


# Optional ABC convenience base classes.
#
# Plugin authors MAY inherit from these to get common boilerplate for free,
# but are not required to do so.
class RdataBase(ABC):
    """Optional ABC for data-reading classes.

    Allows for subclasses to be used as a context manager by providing :meth:`__enter__` 
    and :meth:`__exit__` methods that delegate to :meth:`__del__`. 

    :meth:`__next__` catches :class:`IendError` and converts it to :class:`StopIteration`,
    allowing use of Rdata as an iterator in a for loop or comprehension.

    Subclasses must implement :meth:`__call__`, :meth:`seek_frame`, and :meth:`seek_last`.
    """

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.__del__()

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self.__call__()
        except IendError:
            raise StopIteration

    @abstractmethod
    def __del__(self):
        """Release any resources held by this reader, e.g. open file handles."""

    @abstractmethod
    def __call__(self, nframe: Optional[int] = None):
        """Reads one exposure from the run the :class:`Rdata` is attached to. If
        `nframe` is None, then it will read the frame it is positioned at. If
        nframe is an integer > 0, it will try to read that particular frame;
        if nframe == 0, it reads the last complete frame. nframe == 1 gives
        the first frame. If all works, an MCCD object is returned. It raises
        an HendError if it reaches the end of a local disk file. If it fails
        to read from the server it returns None, but does not raise an
        Exception in order to allow continued attempts as reading from what
        might be a growing file.

        It maintains an attribute 'nframe' corresponding to the frame that the
        file pointer is positioned to read next (most relevant to reading of a
        local file). If it is set to read the last file and that file doesn't
        change it returns None as a sign of failure.

        Parameters:

           nframe : int
              frame number to get, starting at 1. 0 for the last (complete)
              frame. 'None' indicates that the next frame is wanted, unless
              self.nframe = 0 in which case it will try to get the most recent
              frame, whatever that is.

        Apart from reading the raw bytes, the main job of this routine is to
        divide up and re-package the bytes read into Windows suitable for
        constructing CCD objects.

        If access via a server is requested, it is assumed that the file being
        accessed could be being added to.
        """


# ---------------------------------------------------------------------------
# Plugin discovery and loading
# ---------------------------------------------------------------------------
_REQUIRED_ATTRS = ("Rdata", "get_ccd_pars", "supports_server_access")


def discover_instruments() -> dict:
    """Discover all registered instrument plugins.

    Returns a :class:`dict` mapping instrument name to its
    :class:`importlib.metadata.EntryPoint` object.

    Built-in instruments (``hipercam``, ``ultracam``, ``ultraspec``) are
    always present when the hipercam package is installed, because they are
    registered in hipercam's own ``pyproject.toml``.

    Example::

        from hipercam.instruments import discover_instruments
        for name, ep in discover_instruments().items():
            print(name, ep)
    """
    if sys.version_info >= (3, 9):
        from importlib.metadata import entry_points
        eps = entry_points(group="hipercam.instruments")
    else:
        try:
            from importlib.metadata import entry_points
            # Python 3.8: entry_points() returns a plain dict
            eps = entry_points().get("hipercam.instruments", [])
        except ImportError:
            from importlib_metadata import entry_points
            eps = entry_points(group="hipercam.instruments")
    return {ep.name: ep for ep in eps}


def available_instruments() -> list:
    """Return a list of available instrument plugins.

    This is useful to provide a list of possible values for "resource" in 
    the command line scripts like reduce, which expands when plugins are
    added.

    Example::

        from hipercam.instruments import available_instruments
        print(available_instruments())
    """
    # start with the list of built-in instruments
    instruments = ["hs", "hl", "us", "ul", "hf"]
    for instrument_name in discover_instruments().keys():
        # don't duplicate the built-in instruments, which are always present as plugins
        if instrument_name in ['ultracam', 'hipercam', 'ultraspec']:
            continue
        # this can fail if the instrument plugin is broken
        # we should warn instead of crashing
        try:
            instrument = load_instrument(instrument_name)
        except HipercamError as e:
            warnings.warn(f"Failed to load instrument plugin {instrument_name!r}: {e}")
            continue
        instruments.append(f"{instrument_name}:local")
        if instrument.supports_server_access:
            instruments.append(f"{instrument_name}:server")

    return instruments


def load_instrument(name: str):
    """Load and validate a named instrument plugin.

    Loads the module registered under the ``hipercam.instruments``
    entry-point group with the given *name*, checks that it exposes all
    required attributes, and returns it.

    Parameters
    ----------
    name : str
        The instrument name as registered in the entry-point group,
        e.g. ``'hipercam'``, ``'ultracam'``, ``'ultraspec'``.

    Returns
    -------
    module
        The loaded instrument module.

    Raises
    ------
    HipercamError
        If no plugin with that name is registered, or if the loaded module
        is missing required attributes.
    """
    instruments = discover_instruments()
    if name not in instruments:
        available = ", ".join(sorted(instruments)) or "(none)"
        raise HipercamError(
            f"No instrument plugin named {name!r} is registered. "
            f"Available: {available}"
        )

    module = instruments[name].load()

    # Perform an explicit attribute check. runtime_checkable isinstance only
    # verifies method names, not class-level annotations, across all Python
    # versions; the explicit loop here is the reliable guard.
    missing = [attr for attr in _REQUIRED_ATTRS if not hasattr(module, attr)]
    if missing:
        raise HipercamError(
            f"Instrument plugin {name!r} does not satisfy InstrumentModuleProtocol. "
            f"Missing attributes: {missing}"
        )

    return module