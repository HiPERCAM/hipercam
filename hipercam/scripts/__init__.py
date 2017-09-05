"""
command line scripts

All the command line scripts in hipercam are accessed through entry points as
function. This allows them to be included more naturally within the overall
package. Note that the commands 'cadd', 'csub', 'cmul' and 'cdiv' are all
carried out by 'carith'.

The list of functions is a complete list of available scripts. See the
documentation on hipercam.cline for how to supply the arguments on the command
line.
"""

from .carith import carith
from .combine import combine
from .grab import grab
from .hplot import hplot
from .makestuff import makedata, makefield
from .rtplot import rtplot
from .stats import stats

__all__ = [
    'carith', 'combine', 'grab', 'hplot',
    'makedata', 'makefield', 'rtplot', 'stats'
    ]
