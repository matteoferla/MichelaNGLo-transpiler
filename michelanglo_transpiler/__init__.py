#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from warnings import warn
if sys.version_info[0] < 3:
    warn("Sorry man, I told you to use Python 3.")

########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

###############################################################

from .locking_singleton_pymol import GlobalPyMOL

###############################################################

from .transpiler_base_mixin import PyMolTranspiler_base
from .transplier_pse_mixin import PyMolTranspiler_PSE
from .transpiler_mod_mixin import PyMolTranspiler_modifier
from .transpiler_io_mixin import PyMolTranspiler_io


class PyMolTranspiler(PyMolTranspiler_base, PyMolTranspiler_io, PyMolTranspiler_modifier, PyMolTranspiler_PSE):
    """
    The class initialises as a near blank object.
    Historically, it was the transpiler object itself,
    most of the functionality of which now is in ``Transpiler.transpile()``.
    Say ``trans = PyMolTranspiler(job=User.get_username(request)).transpile(file=filename, **settings)``.

    So there are three kinds of bound methods here.

    * Methods with a self-contained pymol session
    * Methods that need pymol started
    * Methods that do not use pymol

    The initialised object does not have a pymol session. ``self.pymol``.
    For list of bound attributes see ``transpiler_base_mixin``.
    """



