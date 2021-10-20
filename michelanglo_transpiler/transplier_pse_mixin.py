########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo  --- PSE Parsing Mixin for PyMolTranspiler
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

########################################################################################################################

from typing import Sequence, Dict, List, Union, Set
import os, re
from copy import deepcopy
from collections import defaultdict
from warnings import warn

import numpy as np
from mako.template import Template
import pymol2


###############################################################


class PyMolTranspiler_PSE:
    ion_names = ('NA', 'MG', 'MN', 'FE', 'CO', '3CO', 'NI', 'NI1', 'NI2', 'CU', 'ZN', 'K', 'SO4', 'PO4', 'CL', 'F', 'BR')
    # The representation of non-bounded spheres appears a bit random

    template_folder = ''  # michelanglo app proper has these.

    def transpile(self, file, view=None, representation=None, skip_disabled=True, combine_objects=True, **settings):
        """
        method that does the conversion of the PSE files.
        Historically view and representation could be set manually. This is likely to not work.
        And is simply a question of requiring the refactoring the code.

        For views see ``.convert_view(view_string)``, which processes the output of PyMOL command `set_view`
        For representation see ``.convert_reps(reps_string)``, which process the output of PyMOL command
        `iterate 1UBQ, print resi, resn,name,ID,reps`

        Combine_objects...
        This is a peculiar/backwards operation to accomodate the case where the objects are not to be merged.
        Whereas normally raw_pdb and ss are str, here they will be dictionaries unfortunately.

        The problem lies with skipping fix_structure() will result in multiple object:molecule

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            # fix structure requires signeton
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            if file:
                self.log(f'[JOB={self.job}] file {file}')
                assert '.pse' in file.lower(), 'Only PSE files accepted.'
                ##orient
                self.pymol.cmd.load(file)
                self.log(f'[JOB={self.job}] File loaded.')
                v = self.pymol.cmd.get_view()
                self.convert_view(v)
                self.log(f'[JOB={self.job}] View converted.')
                if combine_objects:
                    self.fix_structure()
                self.log(f'[JOB={self.job}] Secondary structure fix applied.')
                names_for_mesh_route = []  # this is for a last ditch attempt.
                names_not_mesh = []
                ## General settings
                self.stick_transparency = float(self.pymol.cmd.get('stick_transparency'))
                self.log(f'[JOB={self.job}] 1')
                self.surface_transparency = float(self.pymol.cmd.get('transparency'))
                self.log(f'[JOB={self.job}] 1')
                self.cartoon_transparency = float(self.pymol.cmd.get('cartoon_transparency'))
                self.log(f'[JOB={self.job}] 1')
                self.sphere_transparency = float(self.pymol.cmd.get('sphere_transparency'))
                self.log(f'[JOB={self.job}] 1')
                self.ribbon_transparency = float(self.pymol.cmd.get('ribbon_transparency'))
                self.log(f'[JOB={self.job}] 1')
                self.fov = float(self.pymol.cmd.get("field_of_view"))
                self.log(f'[JOB={self.job}] 1')
                self.fog = float(self.pymol.cmd.get("fog_start")) * 100
                self.log(f'[JOB={self.job}] 1')
                ### sort the pymol objetcs into relevant methods
                self.log(f'[JOB={self.job}] {self.pymol.cmd.get_names()}')
                for obj_name in self.pymol.cmd.get_names():
                    obj = self.pymol.cmd.get_session(obj_name)['names'][0]
                    self.log(f'[JOB={self.job}] {obj_name}')
                    """
                    https://pymolwiki.org/index.php/Get_session
                    0 => name
                    1 => obj or sele?
                    2 => enabled?
                    3 => reps
                    4 => obj_type
                    5 => obj_data... complicated.
                    6 => group_name
                    """
                    # self.pymol.cmd.get_type(obj_name) equivalents in comments
                    if obj[4] == 1:  # object:molecule
                        if obj[1]:
                            continue  # PyMOL selection has no value.
                        if obj[2] == 0:  # PyMOL disabled
                            if skip_disabled:
                                names_not_mesh.append(obj_name)
                            else:
                                raise NotImplementedError()
                        else:  # enabled
                            """ the attr names in get_model differ slightly from the ones iterate gives.
                             as raw pymol output needs to  be an option and
                             the reps variable differs from flags (not even sure they are the same)
                             the get_model way has been depracated, even though iterate seems more barbarous.
                             here is the old code:
                            data = [{'ID': atom.id,
                                     'chain': atom.chain,
                                     'resi': atom.resi, #resi_number is int, but has offset issues?
                                     'resn': atom.resn, #3letter
                                     'name': atom.name,
                                     'elem': atom.symbol,
                                     'reps': atom.flags,
                                     'color': atom.color_code,
                                     'segi': atom.segi}
                                    for atom in pymol.cmd.get_model(obj_name)]
                            """
                            # the object name will be the variable in JS storing the PDB.
                            name = re.sub(r'[^\w_]', '', obj_name)
                            if re.match(r'^\d', name):
                                name = '_'+name
                            self.log(f'[JOB={self.job}] 1')
                            myspace = {'data': []}  # myspace['data'] is the same as self.atoms
                            self.pymol.cmd.iterate(obj_name, self._iterate_cmd, space=myspace)
                            self.log(f'[JOB={self.job}] iterate')
                            self.convert_representation(myspace['data'], **settings)
                            self.log(f'[JOB={self.job}] converted')
                            self.parse_ss(myspace['data']) # resets self.ss each time.
                            # #TODO check if it can be changed to be self.ss = self.parse_ss
                            self.log(f'[JOB={self.job}] 1')
                            self.pdbblocks[name] = '\n'.join(self.ss)+'\n'+self.pymol.cmd.get_pdbstr(obj_name)
                            self.loadfuns[name] = self.get_loadfun_js(tag_wrapped=True, funname=f'load{name}', **settings)
                            names_not_mesh.append(obj_name)
                    elif obj[4] == 2:  # object:map
                        names_for_mesh_route.append(obj_name)
                    elif obj[4] == 3:  # object:mesh
                        names_for_mesh_route.append(obj_name)
                    elif obj[4] == 4:  # 'object:measurement'
                        if obj[2] == 0:  # PyMOL disabled
                            if skip_disabled:
                                names_not_mesh.append(obj_name)
                            else:
                                raise NotImplementedError()
                        else:
                            list_of_coordinates = obj[5][2][0][1]
                            current_distances = []
                            for pi in range(0, len(list_of_coordinates), 6):
                                coord_A = list_of_coordinates[pi:pi + 3]
                                coord_B = list_of_coordinates[pi + 3:pi + 6]
                                current_distances.append({'atom_A': self.get_atom_id_of_coords(coord_A),
                                                          'atom_B': self.get_atom_id_of_coords(coord_B)})
                            self.distances.append({'pairs': current_distances, 'color': obj[5][0][2]})
                    elif obj[4] == 5:  # no idea
                        continue
                    elif obj[4] == 6:  # object:cgo
                        names_for_mesh_route.append(obj_name)
                    elif obj[4] == 7:  # object:surface
                        names_for_mesh_route.append(obj_name)
                    elif obj[4] == 12:  # object:group
                        continue
                self.log(f'[JOB={self.job}] Reps converted.')
                # save and delete all
                self.describe()
                self.pymol.cmd.delete('all')
                if names_for_mesh_route:
                    if 1 == 0:  ##TODO reimplement
                        """
                        This secion has an issue with the alibi transformation.
                        The coordinate vectors need to be moved by the camera movement probably.
                        """
                        objfile = os.path.join(self.temporary_folder, os.path.split(file)[1].replace('.pse', '.obj'))
                        self.pymol.cmd.save(objfile)
                        self.custom_mesh = PyMolTranspiler.convert_mesh(open(objfile, 'r'))
                        os.remove(objfile)
                    else:
                        self.log(f'[JOB={self.job}] WARNING! Conversion of meshes disabled for now.')
            # LEGACY. Likely broken.
            if view:
                self.convert_view(view)
                self.log(f'[JOB={self.job}] View converted.')
            if representation:
                self.convert_representation(representation, **settings)
                self.log(f'[JOB={self.job}] Reps converted.')
            return self

    def get_atom_id_of_coords(self, coord):
        """
        Returns the pymol atom object correspondng to coord. "Needed" for distance.

        :param coord: [x, y, z] vector
        :return: atom

        **PyMOL session**: dependent.
        """
        for on in self.pymol.cmd.get_names(
                enabled_only=1):  # self.pymol.cmd.get_names_of_type('object:molecule') does not handle enabled.
            if self.pymol.cmd.get_type(on) == 'object:molecule':
                o = self.pymol.cmd.get_model(on)
                if o:
                    for atom in o.atom:
                        if all([atom.coord[i] == c for i, c in enumerate(coord)]):
                            return atom
        else:
            return None

    def convert_view(self, view, **settings):
        """
        Converts a Pymol `get_view` output to a NGL M4 matrix.
        If the output is set to string, the string will be a JS command that will require the object stage to exist.
        fog and alpha not implemented.
        self.pymol.cmd.get("field_of_view"))
        self.pymol.cmd.get("fog_start")
        :param view: str or tuple
        :return: np 4x4 matrix or a NGL string

        **PyMOL session**: independent.
        """
        if isinstance(view, str):
            pymolian = np.array([float(i.replace('\\', '').replace(',', '')) for i in view.split() if
                                 i.find('.') > 0])  # isnumber is for ints
        else:
            pymolian = np.array(view)
        self.rotation = pymolian[0:9].reshape([3, 3])
        depth = pymolian[9:12]
        self.z = abs(depth[2]) * 1  # 1 arbitrary weidht to correct. fov should be the same.
        self.position = pymolian[12:15]
        self.teleposition = np.matmul(self.rotation, -depth) + self.position
        self.slab_near = pymolian[11] + pymolian[15]
        self.slab_far = pymolian[11] + pymolian[16]
        # slabs are not clipping...
        self.modrotation = np.multiply(self.rotation, np.array([[-1, -1, -1], [1, 1, 1], [-1, -1, -1]]).transpose())
        c = np.hstack((self.modrotation * self.z, np.zeros((3, 1))))
        m4 = np.vstack((c, np.ones((1, 4))))
        m4[3, 0:3] = -self.position
        self.m4 = m4
        self.validation_text = 'axes\ncgo_arrow [-50,0,0], [50,0,0], gap=0,color=tv_red\n' + \
                               'cgo_arrow [0,-50,0], [0,50,0], gap=0,color=tv_green\n' + \
                               'cgo_arrow [0,0,-50], [0,0,50], gap=0,color=tv_blue\n' + \
                               'cgo_arrow {0}, {1}, gap=0'.format(self.teleposition.tolist(), self.position.tolist()) + \
                               'set_view (\\\n{})'.format(
                                   ',\\\n'.join(['{0:f}, {1:f}, {2:f}'.format(x, y, z) for x, y, z in
                                                 zip(pymolian[:-2:3], pymolian[1:-1:3], pymolian[2::3])]))
        # So it is essential that the numbers be in f format and not e format. or it will be shifted. Likewise for the brackets.
        return self

    def get_view(self, output='matrix', **settings):
        """
        If the output is set to string, the string will be a JS command that will require the object stage to exist.

        :param output: 'matrix' | 'string'
        :return: np 4x4 matrix or a NGL string

        **PyMOL session**: independent.
        """
        warn('This method will be removed soon.', DeprecationWarning)
        assert self.m4 is not None, 'Cannot call get_view without having loaded the data with `convert_view(text)` or loaded a 4x4 transformation matrix (`.m4 =`)'
        if output.lower() == 'string':
            return '//orient\nvar m4 = (new NGL.Matrix4).fromArray({});\nstage.viewerControls.orient(m4);'.format(
                self.m4.reshape(16, ).tolist())
        elif output.lower() == 'matrix':
            return self.m4

    def convert_representation(self, represenation, **settings):
        """iterate all, ID, segi, chain,resi, resn,name, elem,reps, color, ss
        reps seems to be a binary number. controlling the following

        * 0th bit: sticks
        * 7th bit: line
        * 5th bit: cartoon
        * 2th bit: surface

        **PyMOL session**: independent.
        """
        if isinstance(represenation, str):
            text = represenation
            headers = None
            for line in text.split('\n'):
                if not line:
                    continue
                elif line.find('terate') != -1:  # twice. [Ii]terate
                    if line.count(':'):
                        headers = ('ID', 'segi', 'chain', 'resi', 'resn', 'name', 'elem', 'reps', 'color', 'cartoon',
                                   'label')  # gets ignored if iterate> like is present
                    else:
                        headers = [element.rstrip().lstrip() for element in line.split(',')][1:]
                else:
                    # pymol seems to have two alternative outputs.
                    self.atoms.append(dict(zip(headers,
                                               line.replace('(', '').replace(')', '').replace(',', '').replace('\'',
                                                                                                               '').split())))
        else:
            self.atoms = represenation
        # convert reps field. See https://github.com/matteoferla/MichelaNGLo#primitive-equivalence-table for table.
        rep2name = ['sticks', 'spheres', 'surface', 'label', 'non-bounded spheres', 'cartoon', 'ribbon', 'lines',
                    'mesh', 'dots', 'non-bounded', 'cell', 'putty']
        ## determine shape of protein data
        structure = {}
        for atom in self.atoms:
            if atom['chain'] not in structure:
                structure[atom['chain']] = {}
            if atom['resi'] not in structure[atom['chain']]:
                structure[atom['chain']][atom['resi']] = {}
            structure[atom['chain']][atom['resi']][atom['name']] = False
        ## deetermine values for protein
        repdata = {rep2name[i]: deepcopy(structure) for i in (0, 1, 2, 5, 6, 7, 8, 9, 11, 12)}
        for atom in self.atoms:
            reps = [r == '1' for r in reversed("{0:0>12b}".format(int(atom['reps'])))]
            if atom['type'] == 'HETATM':
                if atom['resn'].strip() in self.ion_names:
                    reps[1] = reps[1] or reps[4]  # hetero spheres fix
                reps[7] = reps[7] or reps[10]  # hetero line fix
            # assert atom['chain'], 'The atom has no chain. This ought to be fixed upstream!'
            for i in (0, 1, 2, 5, 6, 7, 8, 9, 11):
                repdata[rep2name[i]][atom['chain']][atom['resi']][atom['name']] = reps[i]
            # label
            if reps[3]:
                self.label['{resi}:{chain}.{name}'.format(**atom)] = atom['label']
        ## deal with putty.
        for atom in self.atoms:
            if atom['name'] == 'CA' and atom['cartoon'] == 7:  # putty override.
                repdata['putty'][atom['chain']][atom['resi']] = {name: True for name in
                                                                 structure[atom['chain']][atom['resi']]}
                # cartoon off!
                repdata['cartoon'][atom['chain']][atom['resi']] = {name: False for name in
                                                                   structure[atom['chain']][atom['resi']]}
        ##convert and collapse
        for i in (0, 1, 2, 5, 6, 7, 8, 9, 11, 12):
            transdata = []
            rep_name = rep2name[i]
            chain_homo_state = True  # are the chains homogeneously represented as rep_name?
            chain_list = []  # a list in case the chains differ
            for chain in repdata[rep_name]:
                resi_homo_state = True  # are the residues homogeneously represented as rep_name?
                resi_list = []  # a list in case the chain is not homogeneous.
                for resi in repdata[rep_name][chain]:
                    if all(repdata[rep_name][chain][resi].values()):
                        resi_list.append(resi)
                    elif any(repdata[rep_name][chain][resi].values()):  # some/all are present
                        if not all(repdata[rep_name][chain][resi].values()):  # some, but not all are present
                            transdata.extend([f'{resi}:{chain}.{name}' for name in repdata[rep_name][chain][resi] if
                                              repdata[rep_name][chain][resi]])
                            resi_homo_state = False
                    else:  # none are.
                        resi_homo_state = False
                if resi_homo_state:  # no residues differ
                    chain_list.append(f':{chain}')
                else:
                    transdata.extend([f'{resi}:{chain}' for resi in self.collapse_list(resi_list)])
                    chain_homo_state = False
            if chain_homo_state:
                transdata.append('*')
            elif len(chain_list):
                transdata.extend(chain_list)
            else:
                pass
            setattr(self, rep_name, transdata)

        # convert color field
        colorset = defaultdict(list)

        # self.swatch[atom['color']]

        def ddictlist():  # a dict of a dict of a list. simple ae?
            return defaultdict(list)

        def tdictlist():  # a dict of a dict of a dict of a list. simple ae?
            return defaultdict(ddictlist)

        carboncolorset = defaultdict(tdictlist)  # chain -> resi -> color_id -> list of atom ids
        colorset = defaultdict(ddictlist)  # element -> color_id -> list of atom ids
        for atom in self.atoms:
            if atom['elem'] == 'C':
                carboncolorset[atom['chain']][atom['resi']][atom['color']].append(atom['ID'])
            else:
                colorset[atom['elem']][atom['color']].append(atom['ID'])
        self.colors = {'carbon': carboncolorset, 'non-carbon': colorset}
        self.convert_color(**settings)
        return self

    @staticmethod
    def collapse_list(l: Sequence) -> List:
        """
        Given a list of residues makes a list of hyphen range string
        """
        icode_polish = lambda resi: int(re.search(r'(\d+)', resi).group(1))
        l = sorted(map(icode_polish, l))
        if len(l) < 2:
            return l
        parts = []
        start = l[0]
        for i in range(1, len(l)):
            fore = l[i - 1]
            aft = l[i]
            if fore + 1 == aft:
                # contiguous
                continue
            else:
                # break
                parts.append(f'{start}-{fore}' if start != fore else str(start))
                start = aft
        parts.append(f'{start}-{aft}' if start != aft else str(start))
        return parts

    def get_reps(self, inner_tabbed=1, stick='sym_licorice', **settings):  # '^'+atom['chain']
        """
        This method is not used.
        """
        warn('This method will be removed soon.', DeprecationWarning)
        assert self.atoms, 'Needs convert_reps first'
        code = ['//representations', 'protein.removeAllRepresentations();']
        if self.colors:
            color_str = 'color: schemeId,'
        else:
            color_str = ''
        if self.lines:
            code.append('var lines = new NGL.Selection( "{0}" );'.format(' or '.join(self.lines)))
            code.append('protein.addRepresentation( "line", {' + color_str + ' sele: lines.string} );')
        if self.sticks:
            code.append('var sticks = new NGL.Selection( "{0}" );'.format(' or '.join(self.sticks)))
            if stick == 'sym_licorice':
                code.append(
                    'protein.addRepresentation( "licorice", {' + color_str + ' sele: sticks.string, multipleBond: "symmetric"} );')
            elif stick == 'licorice':
                code.append('protein.addRepresentation( "licorice", {' + color_str + ' sele: sticks.string} );')
            elif stick == 'hyperball':
                code.append('protein.addRepresentation( "hyperball", {' + color_str + ' sele: sticks.string} );')
            elif stick == 'ball':
                code.append(
                    'protein.addRepresentation( "ball+stick", {' + color_str + ' sele: sticks.string, multipleBond: "symmetric"} );')
        if self.cartoon:
            code.append('var cartoon = new NGL.Selection( "{0}" );'.format(' or '.join(self.cartoon)))
            code.append(
                'protein.addRepresentation( "cartoon", {' + color_str + ' sele: cartoon.string, smoothSheet: true} );')  # capped does not add arrow heads.
        if self.surface:
            code.append('var surf = new NGL.Selection( "{0}" );'.format(' or '.join(self.surface)))
            code.append('protein.addRepresentation( "surface", {' + color_str + ' sele: surf.string} );')
        return code  # self.indent(code, inner_tabbed)

    def convert_color(self, uniform_non_carbon=False, inner_tabbed=1, **settings):
        """
        determine what colors we have.
        ``{'carbon':carboncolorset,'non-carbon': colorset}``
        """
        self.elemental_mapping = {}
        self.catenary_mapping = {}  # pertaining to chains...
        self.residual_mapping = {}
        self.serial_mapping = {}
        # non-carbon
        for elem in self.colors['non-carbon']:  # element -> color_id -> list of atom ids
            if len(self.colors['non-carbon'][elem]) == 1:
                color_id = list(self.colors['non-carbon'][elem].keys())[0]
                self.elemental_mapping[elem] = self.swatch[color_id].hex
            else:
                colors_by_usage = sorted(self.colors['non-carbon'][elem].keys(),
                                         key=lambda c: len(self.colors['non-carbon'][elem][c]), reverse=True)
                self.elemental_mapping[elem] = self.swatch[colors_by_usage[0]].hex
                if not uniform_non_carbon:
                    for color_id in colors_by_usage[1:]:
                        for serial in self.colors['non-carbon'][elem][color_id]:
                            self.serial_mapping[serial] = self.swatch[color_id].hex
        # carbon
        for chain in self.colors['carbon']:
            colors_by_usage = sorted(
                set([col for resi in self.colors['carbon'][chain] for col in self.colors['carbon'][chain][resi]]),
                key=lambda c: len([self.colors['carbon'][chain][resi][c] for resi in self.colors['carbon'][chain] if
                                   c in self.colors['carbon'][chain][resi]]), reverse=True)
            self.catenary_mapping[chain] = self.swatch[colors_by_usage[0]].hex
            for resi in self.colors['carbon'][chain]:  # -> resi -> color_id -> list of atom ids
                if len(self.colors['carbon'][chain][resi]) == 1:
                    color_id = list(self.colors['carbon'][chain][resi].keys())[0]
                    if color_id != colors_by_usage[0]:
                        self.residual_mapping[chain + resi] = self.swatch[color_id].hex
                else:
                    # residue with different colored carbons!
                    for color_id in self.colors['carbon'][chain][resi]:
                        for serial in self.colors['carbon'][chain][resi][color_id]:
                            self.serial_mapping[serial] = self.swatch[color_id].hex
        return self

    def describe(self) -> Dict:
        """
        determine how and what the chains are labelled and what are their ranges.
        ``{'peptide': [f'{first_resi}-{last_resi}:{chain}', ..], 'hetero': [f'[{resn}]{resi}:{chain}', ..]}``

        :rtype: dict

        **PyMOL session**: dependent.
        """
        first_resi = defaultdict(lambda: 9999)
        last_resi = defaultdict(lambda: -9999)
        heteros = set()
        icode_polish = lambda resi: int(re.search(r'(\d+)', resi).group(1))
        for on in self.pymol.cmd.get_names(
                enabled_only=1):  # pymol.cmd.get_names_of_type('object:molecule') does not handle enabled.
            if self.pymol.cmd.get_type(on) == 'object:molecule':
                o = self.pymol.cmd.get_model(on)
                if o:
                    for at in o.atom:
                        if not at.hetatm:
                            if at.resi.isdigit():
                                r = icode_polish(at.resi)
                            else:  ## likely a weird internal residue
                                continue
                            if r < first_resi[at.chain]:
                                first_resi[at.chain] = r
                            if r > last_resi[at.chain]:
                                last_resi[at.chain] = r
                        else:
                            heteros.add((f'{at.resn} and :{at.chain}', None))
        self.description = {
            'peptide': [(f'{first_resi[chain]}-{last_resi[chain]}:{chain}', None) for chain in first_resi],
            'hetero': list(heteros)}
        self.log(f'[JOB={self.job}] description generated.')
        return self.description

    def get_new_letter(self):
        possible = iter('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')
        chains = self.get_chains()
        for c in possible:
            if c not in chains:
                yield c

    def get_chains(self, obj=None) -> Set:
        chains = set()
        if obj == None:
            it = self.pymol.cmd.get_names(enabled_only=1)
        else:
            it = [obj]
        for on in it:
            if self.pymol.cmd.get_type(on) == 'object:molecule':
                o = self.pymol.cmd.get_model(on)
                if o:
                    for atom in o.atom:
                        chains.add(atom.chain)
        return chains

    def get_html(self, ngl='https://cdn.rawgit.com/arose/ngl/v0.10.4-1/dist/ngl.js', **settings):
        """
        Returns a string to be copy-pasted into HTML code.
        :param ngl: (optional) the address to ngl.js. If unspecified it gets it from the RawGit CDN
        :param viewport: (optional) the id of the viewport div, without the hash.
        :param image: (optional) advanced mode with clickable image?
        :return: a string.
        """
        if ngl:
            ngl_string = '<script src="{0}" type="text/javascript"></script>\n'.format(ngl)
        else:
            ngl_string = ''
        if not self.code:
            self.get_js(**settings)
        return '<!-- **inserted code**  -->\n{ngl_string}<script type="text/javascript">{js}</script>\n<!-- **end of code** -->'.format(
            ngl_string=ngl_string,
            js=self.code)

    def write_hmtl(self, template_file='test.mako', output_file='test_generated.html', **kargs):
        if self.verbose:
            print('Making file {0} using template {1}'.format(output_file, template_file))
        template = Template(filename=template_file, format_exceptions=True)
        open(output_file, 'w', newline='\n').write(
            template.render_unicode(transpiler=self, **kargs))
        return self

    def get_js(self, **settings):
        code = Template(filename=os.path.join(self.template_folder, 'output.js.mako')) \
            .render_unicode(structure=self, **settings)
        self.code = code
        return code

    def get_loadfun_js(self, **settings):
        code = Template(filename=os.path.join(self.template_folder, 'loadfun.js.mako')) \
            .render_unicode(structure=self, **settings)
        self.code = code
        return code

    @classmethod
    def convert_mesh(cls, fh, scale=0, centroid_mode='unaltered', origin=None):
        """
        Given a fh or iterable of strings, return a mesh, with optional transformations.
        Note color will be lost.
        Only accepts trianglular meshes!

        :param fh: file handle
        :param scale: 0 do nothing. else Angstrom size
        :param centroid_mode: unaltered | origin | center
        :param origin: if centroid_mode is origin get given a 3d vector.
        :return: {'o_name': object_name, 'triangles': mesh triangles}

        **PyMOL session**: independent. hence why it is a class method.
        """
        mesh = []
        o_name = ''
        scale_factor = 0
        vertices = []
        trilist = []
        sum_centroid = [0, 0, 0]
        min_size = [0, 0, 0]
        max_size = [0, 0, 0]
        centroid = [0, 0, 0]
        for row in fh:
            if row[0] == 'o':
                if o_name:
                    mesh.append({'o_name': o_name, 'triangles': trilist})
                    vertices = []
                    trilist = []
                    scale_factor = 0
                    sum_centroid = [0, 0, 0]
                    min_size = [0, 0, 0]
                    max_size = [0, 0, 0]
                o_name = row.rstrip().replace('o ', '')
            elif row[0] == 'v':
                vertex = [float(e) for e in row.split()[1:]]
                vertices.append(vertex)
                for ax in range(3):
                    sum_centroid[ax] += vertex[ax]
                    min_size[ax] = min(min_size[ax], vertex[ax])
                    max_size[ax] = max(max_size[ax], vertex[ax])
            elif row[0] == 'f':
                if not scale:
                    scale_factor = 1
                elif scale_factor == 0:  # first face.27.7  24.5
                    # euclid = sum([(max_size[ax]-min_size[ax])**2 for ax in range(3)])**0.5
                    scale_factor = scale / max(
                        [abs(max_size[ax] - min_size[ax]) for ax in range(3)])
                    if centroid_mode == 'origin':
                        centroid = [sum_centroid[ax] / len(vertices) for ax in range(3)]
                    elif centroid_mode == 'unaltered':
                        centroid = [0, 0, 0]
                    elif centroid_mode == 'custom':
                        # origin = request.POST['origin'].split(',')
                        centroid = [sum_centroid[ax] / len(vertices) - float(origin[ax]) / scale_factor for ax in
                                    range(3)]  # the user gives scaled origin!
                    else:
                        raise ValueError('Invalid request')
                new_face = [e.split('/')[0] for e in row.split()[1:]]
                if (len(new_face) != 3):
                    pass
                trilist.extend(
                    [int((vertices[int(i) - 1][ax] - centroid[ax]) * scale_factor * 100) / 100 for i in new_face[0:3]
                     for ax in range(3)])
        mesh.append({'o_name': o_name, 'triangles': trilist})
        return mesh
