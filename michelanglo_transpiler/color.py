########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo --- ``GlobalPyMOL`` class
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

########################################################################################################################

from typing import Sequence, Dict, List, Union, Set, Tuple
from pprint import PrettyPrinter

pprint = PrettyPrinter().pprint

###############################################################
from .locking_singleton_pymol import GlobalPyMOL


###############################################################

class ColorItem:
    def __init__(self, value: Sequence):
        """
        ``value`` is a tuple outputed from Pymol by (n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()
        such as ``('bismuth', 5358, (0.6196078658103943, 0.30980393290519714, 0.7098039388656616))``

        :param value: sequence of (name, PyMol index, (R, G, B)) where R, G, B is under 1.
        :type value: sequence of three (str, int, (float, float, float))
        :var name: name of color
        :vartype name: str
        :var index: PyMOL index of color
        :vartype index: int
        :var rgb: R, G, B
        :vartype rgb: Sequence
        :var hex: hex string form of color
        :vartype hex: str
        """
        assert len(
            value) == 3, 'value has to be tuple outputed from Pymol by (n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()'
        self.name = value[0]
        self.index = value[1]
        self.rgb = value[2]
        try:
            self.hex = "0x{0:02x}{1:02x}{2:02x}".format(int(value[2][0] * (2 ** 8 - 1)),
                                                        int(value[2][1] * (2 ** 8 - 1)),
                                                        int(value[2][2] * (2 ** 8 - 1)))
        except TypeError:
            self.hex = "0xffffff"


class ColorSwatch:

    def __init__(self, colors: List[Tuple[str, int, Tuple[float, float, float]]]):
        """
        ColorSwatch()._swatch is a dictionary with indicing being the pymol color number. The values are ColorItem instances.
        Preloading the colors is faster than querying pymol.
        ``print [(n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()]`` in Pymol generates a good amount, but it is not the full amount.

        :param colors: a list like [('white', 0, (1.0, 1.0, 1.0))]
        """
        self._swatch = {}
        for color in colors:
            c = ColorItem(color)
            self._swatch[c.index] = c

    def __getitem__(self, index: int) -> ColorItem:
        """
        :param index: a pymol color index
        """
        if int(index) in self._swatch:
            return self._swatch[int(index)]
        else:
            # warn(f'New color! {index}')
            ci = ColorItem(['', index, GlobalPyMOL.pymol.cmd.get_color_tuple(int(index))])
            self._swatch[ci.index] = ci
            return ci
