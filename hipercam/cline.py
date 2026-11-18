# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Code for handling the command line interface.

This module overrides some of the functions in trm.cline, because that package
is currently not available to make fixes
(see https://github.com/pypi/support/issues/11906).

"""
import signal
import sys

from trm.cline import *


class Cline(Cline):
    """Command line interface for Hipercam scripts."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __exit__(self, *args):
        # This is a temporary fix for a bug in Cline.save().
        # The function tries to prevent SIGINT from interrupting the save,
        # by setting SIGINT to SIG_IGN before the save and restoring it to SIG_DFL
        # after.
        # However, that doesn't actually restore the previous behaviour, which is
        # given by signal.default_int_handler. Notably it means after save() is called
        # SIGINT will no longer raise KeyboardInterrupt, but instead will
        # terminate the process with no exception.
        # That means we could never catch KeyboardInterrupt in the main loop.
        # Since we can't fix trm.cline for now, this is a workaround that should
        # restore the normal behaviour of SIGINT after save() is called.
        # Also, instead of completely suppressing a SIGINT, we'll catch it and
        # raise it again after the save is done. Unlikely to matter since save() should
        # be fast, but it's better practice.
        pending = []
        old_handler = signal.signal(
            signal.SIGINT,
            lambda signum, frame: pending.append(signum)
            )
        try:
            self.save()
        finally:
            signal.signal(
                signal.SIGINT,
                old_handler if old_handler is not None else signal.default_int_handler,
            )
            if pending:
                signal.raise_signal(signal.SIGINT)


    def get_value(self, *args, **kwargs):
        # To make exiting a script nicer while still inputting,
        # catch KeyboardInterrupt and exit cleanly instead of raising an exception.
        try:
            return super().get_value(*args, **kwargs)
        except KeyboardInterrupt:
            print("\nExiting...")
            sys.exit(0)
