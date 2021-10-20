########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo  --- base class for PyMOLTranspiler
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

########################################################################################################################

from typing import Sequence, Dict, List, Union, Set
import os
from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint


from datetime import datetime


###############################################################

from .color import ColorSwatch

###############################################################


class PyMolTranspiler_base:
    # this is not used anymore.
    current_task = f'[{datetime.utcnow()}] idle'
    verbose = False
    swatch = ColorSwatch([('white', 0, (1.0, 1.0, 1.0)), ('black', 1, (0.0, 0.0, 0.0)), ('blue', 2, (0.0, 0.0, 1.0)),
                          ('green', 3, (0.0, 1.0, 0.0)), ('red', 4, (1.0, 0.0, 0.0)),
                          ('cyan', 5, (0.0, 1.0, 1.0)), ('yellow', 6, (1.0, 1.0, 0.0)), ('dash', 7, (1.0, 1.0, 0.0)),
                          ('magenta', 8, (1.0, 0.0, 1.0)),
                          ('salmon', 9, (1.0, 0.6000000238418579, 0.6000000238418579)), ('lime', 10, (0.5, 1.0, 0.5)),
                          ('slate', 11, (0.5, 0.5, 1.0)), ('hotpink', 12, (1.0, 0.0, 0.5)),
                          ('orange', 13, (1.0, 0.5, 0.0)), ('chartreuse', 14, (0.5, 1.0, 0.0)),
                          ('limegreen', 15, (0.0, 1.0, 0.5)), ('purpleblue', 16, (0.5, 0.0, 1.0)),
                          ('marine', 17, (0.0, 0.5, 1.0)),
                          ('olive', 18, (0.7699999809265137, 0.699999988079071, 0.0)),
                          ('purple', 19, (0.75, 0.0, 0.75)), ('teal', 20, (0.0, 0.75, 0.75)),
                          ('ruby', 21, (0.6000000238418579, 0.20000000298023224, 0.20000000298023224)),
                          ('forest', 22, (0.20000000298023224, 0.6000000238418579, 0.20000000298023224)),
                          ('deepblue', 23, (0.25, 0.25, 0.6499999761581421)), ('grey', 24, (0.5, 0.5, 0.5)),
                          ('gray', 25, (0.5, 0.5, 0.5)),
                          ('carbon', 26, (0.20000000298023224, 1.0, 0.20000000298023224)),
                          ('nitrogen', 27, (0.20000000298023224, 0.20000000298023224, 1.0)),
                          ('oxygen', 28, (1.0, 0.30000001192092896, 0.30000001192092896)),
                          ('hydrogen', 29, (0.8999999761581421, 0.8999999761581421, 0.8999999761581421)),
                          ('brightorange', 30, (1.0, 0.699999988079071, 0.20000000298023224)),
                          ('sulfur', 31, (0.8999999761581421, 0.7749999761581421, 0.25)),
                          ('tv_red', 32, (1.0, 0.20000000298023224, 0.20000000298023224)),
                          ('tv_green', 33, (0.20000000298023224, 1.0, 0.20000000298023224)),
                          ('tv_blue', 34, (0.30000001192092896, 0.30000001192092896, 1.0)),
                          ('tv_yellow', 35, (1.0, 1.0, 0.20000000298023224)),
                          ('yelloworange', 36, (1.0, 0.8700000047683716, 0.3700000047683716)),
                          ('tv_orange', 37, (1.0, 0.550000011920929, 0.15000000596046448)),
                          ('pink', 48, (1.0, 0.6499999761581421, 0.8500000238418579)),
                          ('firebrick', 49, (0.6980000138282776, 0.12999999523162842, 0.12999999523162842)),
                          ('chocolate', 50, (0.5550000071525574, 0.22200000286102295, 0.11100000143051147)),
                          ('brown', 51, (0.6499999761581421, 0.3199999928474426, 0.17000000178813934)),
                          ('wheat', 52, (0.9900000095367432, 0.8199999928474426, 0.6499999761581421)),
                          ('violet', 53, (1.0, 0.5, 1.0)),
                          ('lightmagenta', 154, (1.0, 0.20000000298023224, 0.800000011920929)),
                          ('density', 4155, (0.10000000149011612, 0.10000000149011612, 0.6000000238418579)),
                          ('paleyellow', 5256, (1.0, 1.0, 0.5)), ('aquamarine', 5257, (0.5, 1.0, 1.0)),
                          ('deepsalmon', 5258, (1.0, 0.5, 0.5)),
                          ('palegreen', 5259, (0.6499999761581421, 0.8999999761581421, 0.6499999761581421)),
                          ('deepolive', 5260, (0.6000000238418579, 0.6000000238418579, 0.10000000149011612)),
                          ('deeppurple', 5261, (0.6000000238418579, 0.10000000149011612, 0.6000000238418579)),
                          ('deepteal', 5262, (0.10000000149011612, 0.6000000238418579, 0.6000000238418579)),
                          ('lightblue', 5263, (0.75, 0.75, 1.0)), ('lightorange', 5264, (1.0, 0.800000011920929, 0.5)),
                          ('palecyan', 5265, (0.800000011920929, 1.0, 1.0)),
                          ('lightteal', 5266, (0.4000000059604645, 0.699999988079071, 0.699999988079071)),
                          ('splitpea', 5267, (0.5199999809265137, 0.75, 0.0)),
                          ('raspberry', 5268, (0.699999988079071, 0.30000001192092896, 0.4000000059604645)),
                          ('sand', 5269, (0.7200000286102295, 0.550000011920929, 0.30000001192092896)),
                          ('smudge', 5270, (0.550000011920929, 0.699999988079071, 0.4000000059604645)),
                          ('violetpurple', 5271, (0.550000011920929, 0.25, 0.6000000238418579)),
                          ('dirtyviolet', 5272, (0.699999988079071, 0.5, 0.5)),
                          ('deepsalmon', 5273, (1.0, 0.41999998688697815, 0.41999998688697815)),
                          ('lightpink', 5274, (1.0, 0.75, 0.8700000047683716)), ('greencyan', 5275, (0.25, 1.0, 0.75)),
                          ('limon', 5276, (0.75, 1.0, 0.25)),
                          ('skyblue', 5277, (0.20000000298023224, 0.5, 0.800000011920929)),
                          ('bluewhite', 5278, (0.8500000238418579, 0.8500000238418579, 1.0)),
                          ('warmpink', 5279, (0.8500000238418579, 0.20000000298023224, 0.5)),
                          ('darksalmon', 5280, (0.7300000190734863, 0.550000011920929, 0.5199999809265137)),
                          ('helium', 5281, (0.8509804010391235, 1.0, 1.0)),
                          ('lithium', 5282, (0.800000011920929, 0.5019607543945312, 1.0)),
                          ('beryllium', 5283, (0.7607843279838562, 1.0, 0.0)),
                          ('boron', 5284, (1.0, 0.7098039388656616, 0.7098039388656616)),
                          ('fluorine', 5285, (0.7019608020782471, 1.0, 1.0)),
                          ('neon', 5286, (0.7019608020782471, 0.8901960849761963, 0.9607843160629272)),
                          ('sodium', 5287, (0.6705882549285889, 0.3607843220233917, 0.9490196108818054)),
                          ('magnesium', 5288, (0.5411764979362488, 1.0, 0.0)),
                          ('aluminum', 5289, (0.7490196228027344, 0.6509804129600525, 0.6509804129600525)),
                          ('silicon', 5290, (0.9411764740943909, 0.7843137383460999, 0.6274510025978088)),
                          ('phosphorus', 5291, (1.0, 0.5019607543945312, 0.0)),
                          ('chlorine', 5292, (0.12156862765550613, 0.9411764740943909, 0.12156862765550613)),
                          ('argon', 5293, (0.5019607543945312, 0.8196078538894653, 0.8901960849761963)),
                          ('potassium', 5294, (0.5607843399047852, 0.2509803771972656, 0.8313725590705872)),
                          ('calcium', 5295, (0.239215686917305, 1.0, 0.0)),
                          ('scandium', 5296, (0.9019607901573181, 0.9019607901573181, 0.9019607901573181)),
                          ('titanium', 5297, (0.7490196228027344, 0.7607843279838562, 0.7803921699523926)),
                          ('vanadium', 5298, (0.6509804129600525, 0.6509804129600525, 0.6705882549285889)),
                          ('chromium', 5299, (0.5411764979362488, 0.6000000238418579, 0.7803921699523926)),
                          ('manganese', 5300, (0.6117647290229797, 0.47843137383461, 0.7803921699523926)),
                          ('iron', 5301, (0.8784313797950745, 0.4000000059604645, 0.20000000298023224)),
                          ('cobalt', 5302, (0.9411764740943909, 0.5647059082984924, 0.6274510025978088)),
                          ('nickel', 5303, (0.3137255012989044, 0.8156862854957581, 0.3137255012989044)),
                          ('copper', 5304, (0.7843137383460999, 0.5019607543945312, 0.20000000298023224)),
                          ('zinc', 5305, (0.4901960790157318, 0.5019607543945312, 0.6901960968971252)),
                          ('gallium', 5306, (0.7607843279838562, 0.5607843399047852, 0.5607843399047852)),
                          ('germanium', 5307, (0.4000000059604645, 0.5607843399047852, 0.5607843399047852)),
                          ('arsenic', 5308, (0.7411764860153198, 0.5019607543945312, 0.8901960849761963)),
                          ('selenium', 5309, (1.0, 0.6313725709915161, 0.0)),
                          ('bromine', 5310, (0.6509804129600525, 0.16078431904315948, 0.16078431904315948)),
                          ('krypton', 5311, (0.3607843220233917, 0.7215686440467834, 0.8196078538894653)),
                          ('rubidium', 5312, (0.43921568989753723, 0.18039216101169586, 0.6901960968971252)),
                          ('strontium', 5313, (0.0, 1.0, 0.0)), ('yttrium', 5314, (0.5803921818733215, 1.0, 1.0)),
                          ('zirconium', 5315, (0.5803921818733215, 0.8784313797950745, 0.8784313797950745)),
                          ('niobium', 5316, (0.45098039507865906, 0.7607843279838562, 0.7882353067398071)),
                          ('molybdenum', 5317, (0.3294117748737335, 0.7098039388656616, 0.7098039388656616)),
                          ('technetium', 5318, (0.23137255012989044, 0.6196078658103943, 0.6196078658103943)),
                          ('ruthenium', 5319, (0.1411764770746231, 0.5607843399047852, 0.5607843399047852)),
                          ('rhodium', 5320, (0.03921568766236305, 0.4901960790157318, 0.5490196347236633)),
                          ('palladium', 5321, (0.0, 0.4117647111415863, 0.5215686559677124)),
                          ('silver', 5322, (0.7529411911964417, 0.7529411911964417, 0.7529411911964417)),
                          ('cadmium', 5323, (1.0, 0.8509804010391235, 0.5607843399047852)),
                          ('indium', 5324, (0.6509804129600525, 0.4588235318660736, 0.45098039507865906)),
                          ('tin', 5325, (0.4000000059604645, 0.5019607543945312, 0.5019607543945312)),
                          ('antimony', 5326, (0.6196078658103943, 0.38823530077934265, 0.7098039388656616)),
                          ('tellurium', 5327, (0.8313725590705872, 0.47843137383461, 0.0)),
                          ('iodine', 5328, (0.5803921818733215, 0.0, 0.5803921818733215)),
                          ('xenon', 5329, (0.25882354378700256, 0.6196078658103943, 0.6901960968971252)),
                          ('cesium', 5330, (0.34117648005485535, 0.09019608050584793, 0.5607843399047852)),
                          ('barium', 5331, (0.0, 0.7882353067398071, 0.0)),
                          ('lanthanum', 5332, (0.43921568989753723, 0.8313725590705872, 1.0)),
                          ('cerium', 5333, (1.0, 1.0, 0.7803921699523926)),
                          ('praseodymium', 5334, (0.8509804010391235, 1.0, 0.7803921699523926)),
                          ('neodymium', 5335, (0.7803921699523926, 1.0, 0.7803921699523926)),
                          ('promethium', 5336, (0.6392157077789307, 1.0, 0.7803921699523926)),
                          ('samarium', 5337, (0.5607843399047852, 1.0, 0.7803921699523926)),
                          ('europium', 5338, (0.3803921639919281, 1.0, 0.7803921699523926)),
                          ('gadolinium', 5339, (0.2705882489681244, 1.0, 0.7803921699523926)),
                          ('terbium', 5340, (0.1882352977991104, 1.0, 0.7803921699523926)),
                          ('dysprosium', 5341, (0.12156862765550613, 1.0, 0.7803921699523926)),
                          ('holmium', 5342, (0.0, 1.0, 0.6117647290229797)),
                          ('erbium', 5343, (0.0, 0.9019607901573181, 0.4588235318660736)),
                          ('thulium', 5344, (0.0, 0.8313725590705872, 0.32156863808631897)),
                          ('ytterbium', 5345, (0.0, 0.7490196228027344, 0.21960784494876862)),
                          ('lutetium', 5346, (0.0, 0.6705882549285889, 0.1411764770746231)),
                          ('hafnium', 5347, (0.3019607961177826, 0.7607843279838562, 1.0)),
                          ('tantalum', 5348, (0.3019607961177826, 0.6509804129600525, 1.0)),
                          ('tungsten', 5349, (0.12941177189350128, 0.5803921818733215, 0.8392156958580017)),
                          ('rhenium', 5350, (0.14901961386203766, 0.4901960790157318, 0.6705882549285889)),
                          ('osmium', 5351, (0.14901961386203766, 0.4000000059604645, 0.5882353186607361)),
                          ('iridium', 5352, (0.09019608050584793, 0.3294117748737335, 0.529411792755127)),
                          ('platinum', 5353, (0.8156862854957581, 0.8156862854957581, 0.8784313797950745)),
                          ('gold', 5354, (1.0, 0.8196078538894653, 0.13725490868091583)),
                          ('mercury', 5355, (0.7215686440467834, 0.7215686440467834, 0.8156862854957581)),
                          ('thallium', 5356, (0.6509804129600525, 0.3294117748737335, 0.3019607961177826)),
                          ('lead', 5357, (0.34117648005485535, 0.3490196168422699, 0.3803921639919281)),
                          ('bismuth', 5358, (0.6196078658103943, 0.30980393290519714, 0.7098039388656616)),
                          ('polonium', 5359, (0.6705882549285889, 0.3607843220233917, 0.0)),
                          ('astatine', 5360, (0.4588235318660736, 0.30980393290519714, 0.2705882489681244)),
                          ('radon', 5361, (0.25882354378700256, 0.5098039507865906, 0.5882353186607361)),
                          ('francium', 5362, (0.25882354378700256, 0.0, 0.4000000059604645)),
                          ('radium', 5363, (0.0, 0.4901960790157318, 0.0)),
                          ('actinium', 5364, (0.43921568989753723, 0.6705882549285889, 0.9803921580314636)),
                          ('thorium', 5365, (0.0, 0.729411780834198, 1.0)),
                          ('protactinium', 5366, (0.0, 0.6313725709915161, 1.0)),
                          ('uranium', 5367, (0.0, 0.5607843399047852, 1.0)),
                          ('neptunium', 5368, (0.0, 0.5019607543945312, 1.0)),
                          ('plutonium', 5369, (0.0, 0.41960784792900085, 1.0)),
                          ('americium', 5370, (0.3294117748737335, 0.3607843220233917, 0.9490196108818054)),
                          ('curium', 5371, (0.47058823704719543, 0.3607843220233917, 0.8901960849761963)),
                          ('berkelium', 5372, (0.5411764979362488, 0.30980393290519714, 0.8901960849761963)),
                          ('californium', 5373, (0.6313725709915161, 0.21176470816135406, 0.8313725590705872)),
                          ('einsteinium', 5374, (0.7019608020782471, 0.12156862765550613, 0.8313725590705872)),
                          ('fermium', 5375, (0.7019608020782471, 0.12156862765550613, 0.729411780834198)),
                          ('mendelevium', 5376, (0.7019608020782471, 0.05098039284348488, 0.6509804129600525)),
                          ('nobelium', 5377, (0.7411764860153198, 0.05098039284348488, 0.529411792755127)),
                          ('lawrencium', 5378, (0.7803921699523926, 0.0, 0.4000000059604645)),
                          ('rutherfordium', 5379, (0.800000011920929, 0.0, 0.3490196168422699)),
                          ('dubnium', 5380, (0.8196078538894653, 0.0, 0.30980393290519714)),
                          ('seaborgium', 5381, (0.8509804010391235, 0.0, 0.2705882489681244)),
                          ('bohri', 5382, (0.8784313797950745, 0.0, 0.21960784494876862)),
                          ('hassium', 5383, (0.9019607901573181, 0.0, 0.18039216101169586)),
                          ('meitnerium', 5384, (0.9215686321258545, 0.0, 0.14901961386203766)),
                          ('deuterium', 5385, (0.8999999761581421, 0.8999999761581421, 0.8999999761581421)),
                          ('lonepair', 5386, (0.5, 0.5, 0.5)),
                          ('pseudoatom', 5387, (0.8999999761581421, 0.8999999761581421, 0.8999999761581421))])
    temporary_folder = 'temp' # overridden in michelanglo_app.setup_folders
    _iterate_cmd = "data.append({'ID': ID,  'segi': segi, 'chain': chain, 'resi': resi, 'resn': resn, 'name':name, 'elem':elem, 'reps':reps, 'color':color, 'ss': ss, 'cartoon': cartoon, 'label': label, 'type': type})"
    boring_ligand = (  # 'WAT', 'HOH',  # `TP3` water is ambiguous and rare
        'LI', 'NA', 'K', 'RB',  # group 1 cations
        'BE', 'MG', 'CA', 'SR',  # earth metal cations
        'F', 'CL', 'BR', 'I',  # halogens
        'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',  # top period transition metals
        '3CO',  # cobalt (iii) ion
        'BUQ',  # 4-hydroxy-2-butanone
        # 'NAG',  # n-acetyl-d-glucosamine
        # 'NAD',  # nicotinamide-adenine-dinucleotide
        'CR',  # chromium ion
        # 'SF4',  # iron/sulfur cluster
        'EOH',  # ethanol
        'ZNO',  # zinc ion, 2 waters coordinated
        'NAO',  # sodium ion, 1 water coordinated
        'EOM',  # ethyloxymethoxyl
        'EHN',  # ethane
        # 'NAP',  # nadp nicotinamide-adenine-dinucleotide phosphate
        'CCN',  # acetonitrile
        'NAW',  # sodium ion, 3 waters coordinated
        'BR',  # bromide ion
        'EGL',  # ethylene glycol
        'NI2',  # nickel (ii) ion, 2 waters coordinated
        # 'GSH',  # glutathione
        'NI1',  # nickel ion, 1 water coordinated
        # 'O2',  # oxygen molecule
        'BA',  # barium ion
        'RU',  # ruthenium ion
        # 'SAH',  # s-adenosyl-l-homocysteine
        'GU7',  # 2-amino-7-[2-(2-hydroxy-1-hydroxymethyl-ethylamino)-ethyl]-1,7-dihydro-purin-6-one
        # 'SAM',  # s-adenosylmethionine
        'TAS',  # trihydroxyarsenite(iii)
        'DCE',  # 1,2-dichloroethane
        '2BM',  # dibromomethane
        # 'TP7',  # coenzyme b
        'OF3',  # ferric ion, 1 water coordinated
        'OF1',  # ferrous ion, 1 water coordinated
        'RB',  # rubidium ion
        'IOH',  # 2-propanol, isopropanol
        'MW1',  # manganese ion, 1 water coordinated
        'IOD',  # iodide ion
        'C2O',  # cu-o-cu linkage
        'BNZ',  # benzene
        'TCN',  # tetracyanonickelate ion
        'ARS',  # arsenic
        'NH4',  # ammonium ion
        'GD',  # gadolinium atom
        # 'PER',  # peroxide ion
        'GA',  # gallium (iii) ion
        # 'TPP',  # thiamine diphosphate
        'CHX',  # cyclohexane
        # 'CME',  # s,s-(2-hydroxyethyl)thiocysteine
        # 'THB',  # tetrahydrobiopterin
        'IPA',  # isopropyl alcohol
        'CD1',  # cadmium ion, 1 water coordinated
        'OH',  # hydroxide ion
        'SO4',  # sulfate ion
        'DTT',  # 2,3-dihydroxy-1,4-dithiobutane
        # 'PQN',  # phylloquinone
        'CYN',  # cyanide ion
        # 'PQQ',  # pyrroloquinoline quinone
        'PYJ',  # phenylethane
        # 'PEO',  # hydrogen peroxide
        'NA6',  # sodium ion, 6 waters coordinated
        'MBR',  # tribromomethane
        'NA5',  # sodium ion, 5 waters coordinated
        'OS',  # osmium ion
        'MAN',  # alpha-d-mannose
        'CMO',  # carbon monoxide
        'OCL',  # cobalt ion, 1 water coordinated
        'DMF',  # dimethylformamide
        'OCN',  # cobalt ion, 2 waters coordinated
        'MO3',  # magnesium ion, 3 waters coordinated
        'NGN',  # nitrogen
        'ACT',  # acetate ion
        'U1',  # uranium atom
        'HDZ',  # nitrogen molecule
        'MO5',  # magnesium ion, 5 waters coordinated
        'MO4',  # magnesium ion, 4 waters coordinated
        'VO4',  # vanadate ion
        'DMS',  # dimethyl sulfoxide
        'FUC',  # alpha-l-fucose
        'PCL',  # platinum(ii) di-chloride
        'CB5',  # cobalt bis(1,2-dicarbollide)
        'EEE',  # ethyl acetate
        'HG',  # mercury (ii) ion
        'NO2',  # nitrite ion
        # 'CMP',  # adenosine-3',5'-cyclic-monophosphate
        'PR',  # praseodymium ion
        'BMA',  # beta-d-mannose
        'IUM',  # uranyl (vi) ion
        'PT',  # platinum (ii) ion
        'ZN2',  # zinc ion on 3-fold crystal axis
        # 'TTP',  # thymidine-5'-triphosphate
        'NO3',  # nitrate ion
        'YT3',  # yttrium (iii) ion
        # 'TYS',  # o-sulfo-l-tyrosine
        'PB',  # lead (ii) ion
        'M2M',  # 1-methoxy-2-(2-methoxyethoxy)ethane
        'ZO3',  # zinc ion, 3 waters coordinated
        'PD',  # palladium ion
        # 'AMP',  # adenosine monophosphate
        'PI',  # hydrogenphosphate ion
        'MH3',  # manganese ion, 1 hydroxyl coordinated
        'AF3',  # aluminum fluoride
        'ZN',  # zinc ion
        'MN3',  # manganese (iii) ion
        'OXY',  # oxygen molecule
        'NI',  # nickel (ii) ion
        # 'CSD',  # 3-sulfinoalanine
        # 'OX',  # bound oxygen
        'PS5',  # pentasulfide-sulfur
        'MN5',  # manganese ion, 5 waters coordinated
        'MN6',  # manganese ion, 6 waters coordinated
        'S',  # sulfur atom
        'HOH',  # water
        'W',  # tungsten ion
        'SB',  # antimony (iii) ion
        # 'FOL',  # folic acid
        'OXE',  # ortho-xylene
        'PT4',  # platinum (iv) ion
        'PBM',  # trimethyl lead ion
        'O',  # oxygen atom
        'MW2',  # manganese dihydrate ion
        'MG',  # magnesium ion
        '543',  # calcium ion, 6 waters plus ethanol coordinated
        'MSM',  # (methylsulfanyl)methane
        # 'C5P',  # cytidine-5'-monophosphate
        'ANL',  # aniline
        'MTO',  # bound water
        'NO',  # nitric oxide
        'TBU',  # tertiary-butyl alcohol
        'OPY',  # (3s)-4-oxo-4-piperidin-1-ylbutane-1,3-diamine
        'PC4',  # tetrachloroplatinate(ii)
        # 'GU3',  # methyl 3-o-methyl-2,6-di-o-sulfo-alpha-d-glucopyranoside
        # 'GU2',  # 2,3-di-o-methyl-alpha-l-idopyranuronic acid
        # 'GU1',  # 2,3-di-o-methyl-beta-d-glucopyranuronic acid
        'MOH',  # methanol
        # 'ANP',  # phosphoaminophosphonic acid-adenylate ester
        # 'GU6',  # 2,3,6-tri-o-sulfonato-alpha-d-glucopyranose
        # 'GU5',  # 2,3-di-o-methyl-6-o-sulfonato-alpha-d-glucopyranose
        # 'GU4',  # 2,3,4,6-tetra-o-sulfonato-alpha-d-glucopyranose
        'AU',  # gold ion
        'OC3',  # calcium ion, 3 waters coordinated
        'BTN',  # biotin
        'I42',  # hydroxy(dioxido)oxovanadium
        'OC4',  # calcium ion, 4 waters coordinated
        'OC7',  # calcium ion, 7 waters coordinated
        'OC6',  # calcium ion, 6 waters coordinated
        # 'TMP',  # thymidine-5'-phosphate
        'RE',  # rhenium
        'GD3',  # gadolinium ion
        # 'CTP',  # cytidine-5'-triphosphate
        'ACE',  # acetyl group
        '3OF',  # hydrated fe (iii) ion, 2 waters coordinated
        'ETZ',  # diethyl ether
        'MM4',  # molybdenum (iv) oxide
        'IN',  # indium (iii) ion
        'ACN',  # acetone
        'DOD',  # deuterated water
        'AST',  # arsenite
        # 'COA',  # coenzyme a
        'EU',  # europium ion
        'DOX',  # dioxane
        # 'COB',  # co-methylcobalamin
        # 'B12',  # cobalamin
        'REO',  # perrhenate
        # 'ATP',  # adenosine-5'-triphosphate
        'CD3',  # cadmium ion, 3 waters coordinated
        # 'U10',  # ubiquinone-10
        'ACY',  # acetic acid
        'PEG',  # di(hydroxyethyl)ether
        'YB',  # ytterbium (iii) ion
        # 'NDP',  # nadph dihydro-nicotinamide-adenine-dinucleotide phosphate
        'NBZ',  # nitrobenzene
        'ETI',  # iodoethane
        'C2C',  # cu-cl-cu linkage
        'NA',  # sodium ion
        'FMT',  # formic acid
        'ASC',  # ascorbic acid
        'AU3',  # gold 3+ ion
        'FE2',  # fe (ii) ion
        'LNK',  # pentane
        'SEK',  # selenocyanate ion
        'MO1',  # magnesium ion, 1 water coordinated
        'EU3',  # europium (iii) ion
        '1BO',  # 1-butanol
        'AUC',  # gold (i) cyanide ion
        'CLO',  # chloro group
        'FE',  # fe (iii) ion
        'DUM',  # dummy atoms
        # 'ADP',  # adenosine-5'-diphosphate
        'OF2',  # 2 ferric ion, 1 bridging oxygen
        'BEF',  # beryllium trifluoride ion
        'FEL',  # hydrated fe
        'BF4',  # beryllium tetrafluoride ion
        'HEX',  # hexane
        'CUZ',  # (mu-4-sulfido)-tetra-nuclear copper ion
        # 'NDG',  # 2-(acetylamino)-2-deoxy-a-d-glucopyranose
        'XE',  # xenon
        # 'FMN',  # flavin mononucleotide
        'YAN',  # 1,2-dichlorobenzene
        'CUA',  # dinuclear copper ion
        'V',  # vanadium ion
        'CUO',  # cu2-o2 cluster
        # 'HEM',  # protoporphyrin ix containing fe
        # 'GMP',  # guanosine
        'CU',  # copper (ii) ion
        'MGF',  # trifluoromagnesate
        # 'GDP',  # guanosine-5'-diphosphate
        'CFT',  # trifluoromethane
        'SBT',  # 2-butanol
        # 'PLP',  # pyridoxal-5'-phosphate
        'SR',  # strontium ion
        'FU1',  # tetrahydrofuran
        'EDN',  # ethane-1,2-diamine
        'EDO',  # 1,2-ethanediol
        'H2S',  # hydrosulfuric acid
        'ND4',  # ammonium cation with d
        'BRO',  # bromo group
        'KR',  # krypton
        'CS',  # cesium ion
        'NME',  # methylamine
        # 'CDP',  # cytidine-5'-diphosphate
        'HGI',  # mercury (ii) iodide
        'SM',  # samarium (iii) ion
        # 'ALY',  # n(6)-acetyllysine
        # 'NMO',  # nitrogen monoxide
        # 'TDP',  # thiamin diphosphate
        'SE',  # selenium atom
        'HO',  # holmium atom
        '3CN',  # 3-aminopropane
        'AZI',  # azide ion
        # 'F42',  # coenzyme f420
        'FLO',  # fluoro group
        '6MO',  # molybdenum(vi) ion
        'EMC',  # ethyl mercury ion
        'Y1',  # yttrium ion
        # 'MO7', # bis(mu4-oxo)-bis(mu3-oxo)-octakis(mu2-oxo)-dodecaoxo-heptamolybdenum (vi)
        'SE4',  # selenate ion
        'BF2',  # beryllium difluoride
        'CO',  # cobalt (ii) ion
        # 'NGD', # 3-(aminocarbonyl)-1-[(2r,3r,4s,5r)-5-({[(s)-{[(s)-{[(2r,3s,4r,5r)-5-(2-amino-6-oxo-1,6-dihydro-9h-purin-9-yl)-3,4-dihydroxytetrahydrofuran-2-yl]methoxy}(hydroxy)phosphoryl]oxy}(hydroxy)phosphoryl]oxy}methyl)-3,4-dihydroxytetrahydrofuran-2-yl]pyridinium
        '2MO',  # molybdenum (iv)oxide
        '202',  # bromic acid
        'DIS',  # disordered solvent
        'MBN',  # toluene
        'LA',  # lanthanum (iii) ion
        'PGO',  # s-1,2-propanediol
        'CL',  # chloride ion
        'HP6',  # heptane
        'SO2',  # sulfur dioxide
        'LI',  # lithium ion
        # 'PPS',  # 3'-phosphate-adenosine-5'-phosphate sulfate
        # 'TPO',  # phosphothreonine
        'POL',  # n-propanol
        # 'GU0',  # 2,3,6-tri-o-sulfonato-alpha-l-galactopyranose
        'SGM',  # monothioglycerol
        'DTU',  # (2r,3s)-1,4-dimercaptobutane-2,3-diol
        'MOO',  # molybdate ion
        'TE',  # tellurium
        'TB',  # terbium(iii) ion
        'CA',  # calcium ion
        # 'FAD',  # flavin-adenine dinucleotide
        'CNV',  # propanenitrile
        'GOL',  # glycerol
        'SCN',  # thiocyanate ion
        'AG',  # silver ion
        'PO4',  # phosphate ion
        'IR',  # iridium ion
        'DIO',  # 1,4-diethylene dioxide
        'NH2',  # amino group
        '8CL',  # chlorobenzene
        '3NI',  # nickel (iii) ion
        'IRI',  # iridium hexammine ion
        # 'UTP',  # uridine 5'-triphosphate
        'AR',  # argon
        # 'N4M', # 5-formyltetrahydromethanopterin
        'CE',  # cerium (iii) ion
        'NH3',  # ammonia
        'MN',  # manganese (ii) ion
        'CNN',  # cyanamide
        'HGC',  # methyl mercury ion
        # 'GU8',  # 2,3,6-tri-o-methyl-beta-d-glucopyranose
        # 'GTP',  # guanosine-5'-triphosphate
        # 'UDP',  # uridine-5'-diphosphate
        'OC2',  # calcium ion, 2 waters coordinated
        'ART',  # arsenate
        'TFH',  # nitrogen of trifluoro-ethylhydrazine
        'MCH',  # trichloromethane
        '2NO',  # nitrogen dioxide
        '6WO',  # oxo-tungsten(vi)
        'CD5',  # cadmium ion, 5 waters coordinated
        # 'KCX',  # lysine nz-carboxylic acid
        'E1H',  # ethanimine
        'ARF',  # formamide
        'TL',  # thallium (i) ion
        'DXE',  # 1,2-dimethoxyethane
        # 'GU9',  # 2,3,6-tri-o-methyl-alpha-d-glucopyranose
        'IDO',  # iodo group
        'KO4',  # potassium ion, 4 waters coordinated
        'NRU',  # ruthenium (iii) hexaamine ion
        '4MO'  # molybdenum(iv) ion
    )
    aa_ligand = (
        'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
        'THR', 'VAL', 'TRP', 'TYR')

    water_ligand = ('HOH', 'WAT', 'TP3')

    def __init__(self, verbose: bool = False,
                 validation: bool = False,
                 pdb: str = '',
                 job: str = 'task'):
        """
        Converter. ``__init__`` does not interact with PyMOL.

        :param: job: this is needed for the async querying of progress in the app, but not the transpiler code itself. see .log method
        :param: file: filename of PSE file.
        :param verbose: print?
        :param validation: print validation_text set for pymol?
        :param view: the text from PymOL get_view
        :param representation: the text from PyMOL iterate
        :param pdb: the PDB name or code
        """
        self.job = job
        self.verbose = verbose
        self.validation = validation  # boolean for printing.
        self.validation_text = ''
        self.pdb = pdb # False(/None?) or 4 letter code.
        self.pdbblocks = {} # the pymol transpiler now uses this, but everything else does not.
        self.loadfuns = {}
        self.rotation = None
        self.modrotation = None
        self.position = None
        self.teleposition = None
        self.scale = 10
        self.slab_far = 1
        self.slab_near = 1
        self.fov = 40
        self.m4 = None
        self.m4_alt = None
        self.notes = ''
        self.atoms = []
        self.sticks = []  # 0 licorice
        self.stick_transparency = 0
        self.spheres = []  # 1 spacefill
        self.sphere_transparency = 0
        self.surface = []  # 2 surface
        self.surface_transparency = 0
        self.label = {}  # sele : text
        self.cartoon = []  # 5 cartoon
        self.cartoon_transparency = 0
        self.ribbon = []  # 6 backbone
        self.ribbon_transparency = 0
        self.lines = []  # 7 line
        self.mesh = []  # 8 surface {contour: true}
        self.dots = []  # 9 point
        self.cell = []  # 11 cell
        self.putty = []  # NA.
        self.colors = []
        self.distances = []  # and h-bonds...
        self.ss = []
        self.headers = []  # overwrides .ss in dynamic attr pdb_block
        self.code = ''
        self.pymol = None
        self.raw_pdb = None  # this is set from the instance `prot.raw_pdb = open(file).read()`
        self.custom_mesh = []
        self.description = {}

    @classmethod
    def log(cls, msg: str):
        """
        Lazy logging. Only logs the last event.
        """
        if cls.verbose:
            print(f'DEBUG {msg}')
        cls.current_task = f'[{datetime.utcnow()} GMT] {msg}'
