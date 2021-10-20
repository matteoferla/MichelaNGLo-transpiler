########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo  --- Mixins that mods pdbs for PyMolTranspiler
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

########################################################################################################################

from typing import *
import os, re
import pymol2
from Bio.Data.IUPACData import protein_letters_1to3 as p1to3
import logging


log = logging.getLogger()


###############################################################

from .locking_singleton_pymol import GlobalPyMOL

###############################################################

class PyMolTranspiler_modifier:

    def renumber(self,
                 pdb:str,
                 definitions:List,
                 make_A:Union[str,None]=None,
                 remove_solvent:bool=True,
                 sequence: Optional[str]=None):
        """
        Fetches a pdb file into a transpiler obj.
        The sequence was an attempt to fix a bug that was a mistake in a test. it is not needed for production.

        :param file: str file name
        :param definitions: Structure.chain_definitions e.g. [{'chain': 'A', 'uniprot': 'Q9BZ29', 'x': 1605, 'y': 2069, 'offset': 1604, 'range': '1605-2069', 'name': None, 'description': None},
        :return: self

        **PyMOL session**: self-contained.
        """
        log.debug('Renumbering...')
        assert pdb is not None, 'No PDB block provided'
        assert pdb != '', 'Blank PDB block provided'
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            if len(pdb) == 4: ##Is this block redundant?
                self.pymol.cmd.fetch(pdb, type='pdb')  ## using PDB for simplicity. Using CIF may be nicer...
            else:
                self.pymol.cmd.read_pdbstr(pdb,'blockprotein')
            for chain in definitions:
                if chain["offset"] != 0:
                    #print(f'chain {chain["chain"]}', f'resi=str(int(resi){chain["offset"]:+d})')
                    self.pymol.cmd.alter(f'chain {chain["chain"]}', f'resv += {chain["offset"]}')
                    self.pymol.cmd.sort()
            if remove_solvent:
                self.pymol.cmd.remove('solvent')
            if make_A is not None and make_A != 'A':
                log.debug(self.pymol.cmd.get_fastastr('chain A'))
                self.pymol.cmd.alter('chain A', 'chain="XXX"')
                self.pymol.cmd.sort()
                log.debug(self.pymol.cmd.get_fastastr('chain XXX'))
                self.pymol.cmd.alter(f'chain {make_A}', 'chain="A"')
                self.pymol.cmd.sort()
                log.debug(self.pymol.cmd.get_fastastr('chain A'))
                self.pymol.cmd.alter('chain XXX', f'chain="{make_A}"')
                log.debug(self.pymol.cmd.get_fastastr(f'chain {make_A}'))
                self.pymol.cmd.sort()
            self.fix_structure()
            self.pymol.cmd.sort()
            if sequence:
                # assumes
                # * chain A
                # * it is a shifted issue due to missing residue making a shifted offset
                # * that there is not a residue repeat
                # make_A is chain letter in the definitions but A in the structure
                chain_A = [chain for chain in definitions if chain['chain'] == make_A][0]

                x = chain_A['x']
                if self.pymol.cmd.select('chain A') == 0:
                    raise ValueError('Somehow lost chain A.')
                pdb_seq = self.pymol.cmd.get_fastastr(f'resi {x}-{x+19} and chain A').split('\n')[1]
                ref_seq = sequence[x - 1:x+20]
                if (pdb_seq[0] != ref_seq[0]) or (pdb_seq[1] != ref_seq[1]):
                    for shift in range(20):
                        if (pdb_seq[0] == ref_seq[0+shift]) and ((pdb_seq[1] == ref_seq[1+shift])):
                            # self.pymol.cmd.alter('chain A', f'resv+={shift}')
                            # self.pymol.cmd.sort()
                            # # debug:
                            # print('shifted', shift)
                            # pdb_seq = self.pymol.cmd.get_fastastr(f'resi {x}-{x + 20} and chain A').split('\n')[1]
                            # assert (pdb_seq[0] != ref_seq[0]) or (pdb_seq[1] != ref_seq[1]), 'Did not shift'
                            break
            self.parse_ss()
            self.raw_pdb = self.remove_anisou(self.pymol.cmd.get_pdbstr())
            return self

    def sdf_to_pdb(self, infile: str, reffile: str) -> str:
        """
        A special class method to convert a sdf to pdb but with the atom index shifted so that the pdb can be cat'ed.

        :param infile: sdf file
        :param reffile: pdb file for the indices.
        :return: PDB block

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            combofile = infile.replace('.sdf', '_combo.pdb')
            minusfile = infile.replace('.sdf', '_ref.pdb')
            self.pymol.cmd.load(infile, 'ligand')
            self.pymol.cmd.alter('all', 'chain="Z"')
            self.pymol.cmd.load(reffile, 'apo')
            self.pymol.cmd.alter('all','segi=""')
            self.pymol.cmd.sort()
            self.pymol.cmd.create('combo','apo or ligand')
            self.pymol.cmd.save(combofile, 'combo')
            self.pymol.cmd.save(minusfile, 'apo')
            with open(minusfile) as fh:
                ref = fh.readlines()
            with open(combofile) as fh:
                combo = fh.readlines()
            ligand = ''.join([line for line in combo if line not in ref and line.strip() != ''])
            os.remove(combofile)
            os.remove(minusfile)
            return ligand

    def _mutagen(self, mutations:List, chain: str = None, chains: List = None) -> None:
        """
        Create a mutant protein based on a list of mutations on the already loaded protein.
        To use the pymol2 module it uses https://github.com/schrodinger/pymol-open-source/issues/76

        :param outfile: str the file to save the mod as.
        :param mutations: list of string in the single letter format (A234P) without "p.".
        :param chain: str chain id in the pdb loaded.
        :param chains: it must be same len as mutation.
        :return: None

        **PyMOL session**: dependent.
        """
        self.pymol.cmd.wizard("mutagenesis", _self=self.pymol.cmd)
        self.pymol.cmd.do("refresh_wizard")
        for i, mutant in enumerate(mutations):
            if chains:
                chain = chains[i]
            mutant = mutant.replace('p.','').strip()
            if re.search("(\d+)", mutant) is None:
                raise ValueError(f'{mutant} is not a mutation - like A24S or p.Ala24Ser')
            n = re.search("(\d+)", mutant).group(1)
            if re.match("\w{3}\d+\w{3}", mutant):  # 3 letter Arg
                f = re.match("\w{3}\d+(\w{3})", mutant).group(1).upper()
            elif re.match("\w{1}\d+\w{1}", mutant) and mutant[-1] in p1to3:  # 1 letter R
                f = p1to3[mutant[-1]].upper()
            else:
                raise ValueError(f'{mutant} is not a valid mutation. It should be like A123W. Truncations not allowed.')
            #print('f looks like ',f)
            #print('sele ',f"{chain}/{n}/")
            self.pymol.cmd.get_wizard().set_mode(f)
            try:
                self.pymol.cmd.get_wizard().do_select(f"{chain}/{n}/")
            except self.pymol.parsing.QuietException: #color. cf. https://github.com/schrodinger/pymol-open-source/issues/76
                pass
            self.pymol.cmd.get_wizard().apply()
            #m = self.pymol.cmd.get_model(f"resi {n} and name CA").atom
            #if m:
            #    pass
            #    # assert f == m[0].resn, f'Something is not right {r} has a {m[0].atom}'

    def mutate_block(self, block:str, mutations: List, chain=None, chains=None):
        """
        Create a mutant protein based on a list of mutations on a PDB code.
        :param block: str pdb block
        :param outfile: str the file to save the mod as.
        :param mutations: list of string in the single letter format (A234P) without "p.".
        :param chain: str chain id in the pdb loaded.
        :return:


        **PyMOL session**: self-contained.
        """
        with GlobalPyMOL() as self.pymol:
            self.pymol.cmd.read_pdbstr(block,'blockprotein')
            self.headers = [row.replace('"','').replace("'",'').replace("\\",'') for row in block.split("\n") if any([k in row for k in ('LINK', 'HELIX', 'SHEET')])]
            self._mutagen(mutations=mutations, chain=chain, chains=chains)
            self.raw_pdb = self.pymol.cmd.get_pdbstr()
        return self

    def chain_removal_block(self, block: str, chains : List):
        """
        ** PyMOL session **: self - contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.read_pdbstr(block,'blockprotein')
            for chain in chains:
                self.pymol.cmd.remove(f'chain {chain}')
            self.raw_pdb = self.pymol.cmd.get_pdbstr()
        return self

    def dehydrate_block(self, block: str, water:bool=False, ligand:bool=False):
        """
        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.read_pdbstr(block,'blockprotein')
            if water:
                self.pymol.cmd.remove('solvent')
            if ligand:
                self.pymol.cmd.remove(' or '.join([f'resn {l}' for l in self.boring_ligand]))
            self.raw_pdb = self.pymol.cmd.get_pdbstr()
        return self

    #### To be removed

    def mutate_code(self, code, outfile, mutations, chain=None, chains=None):
        """
        Create a mutant protein based on a list of mutations on a PDB code.
        :param code: str pdb code.
        :param outfile: str the file to save the mod as.
        :param mutations: list of string in the single letter format (A234P) without "p.".
        :param chain: str chain id in the pdb loaded.
        :return:


        **PyMOL session**: self-contained.
        """
        with GlobalPyMOL() as self.pymol:
            self.pymol.cmd.fetch(code)
            self._mutagen(mutations=mutations, chain=chain, chains=chains)
            self.pymol.cmd.save(outfile)
            self.pymol.cmd.delete('all')
        return 1

    def mutate_file(self, infile:str, outfile:str, mutations:List[str], chain:str=None, chains:List=None):
        """
        Create a mutant protein based on a list of mutations on a PDB file path.

        :param infile: str
        :param outfile: str the file to save the mod as.
        :param mutations: list of string in the single letter format (A234P) without "p.".
        :param chain: str chain id in the pdb loaded.
        :return:

        **PyMOL session**: self-contained.
        """
        with GlobalPyMOL() as self.pymol:
            self.pymol.cmd.load(infile)
            self._mutagen(mutations=mutations, chain=chain, chains=chains)
            self.pymol.cmd.save(outfile)
            self.pymol.cmd.delete('all')
        return 1

    def dehydrate_code(self, code:str, outfile:str, water=False, ligand=False):
        """
        Create a mutant protein based on a list of mutations on a PDB code.
        :param code: str pdb code.
        :param outfile: str the file to save the mod as.
        :param mutations: list of string in the single letter format (A234P) without "p.".
        :param chain: str chain id in the pdb loaded.
        :return:

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            self.pymol.cmd.fetch(code)
            if water:
                self.pymol.cmd.remove('solvent')
            if ligand:
                self.pymol.cmd.remove(' or '.join([f'resn {l}' for l in self.boring_ligand]))
            self.pymol.cmd.save(outfile)
            self.pymol.cmd.delete('all')
        return 1

    def dehydrate_file(self, infile:str, outfile:str, water=False, ligand=False):
        """
        Create a mutant protein based on a list of mutations on a PDB file path.

        :param infile: str
        :param outfile: str the file to save the mod as.
        :param mutations: list of string in the single letter format (A234P) without "p.".
        :param chain: str chain id in the pdb loaded.
        :return:

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            self.pymol.cmd.load(infile)
            if water:
                self.pymol.cmd.remove('solvent')
            if ligand:
                self.pymol.cmd.remove(' or '.join([f'resn {l}' for l in self.boring_ligand]))
            self.pymol.cmd.save(outfile)
            self.pymol.cmd.delete('all')
        return 1

    def _chain_removal(self, outfile, chains):
        """
        Create a mutant protein based on a list of mutations on the already loaded protein.

        :param outfile: str the file to save the mod as.
        :param chains: list str chain id in the pdb loaded.
        :return: None

        **PyMOL session**: dependent
        """
        for chain in chains:
            self.pymol.cmd.remove(f'chain {chain}')
        self.pymol.cmd.save(outfile)
        self.pymol.cmd.delete('all')

    def chain_removal_code(self, code, outfile, chains):
        """
        Create a mutant protein based on a list of mutations on a PDB code.
        :param code: str pdb code.
        :param outfile: str the file to save the mod as.
        :param chains: list of str chain id in the pdb loaded.
        :return:

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            self.pymol.cmd.fetch(code)
            self._chain_removal(outfile, chains)
        return 1

    def chain_removal_file(self, infile, outfile, chains):
        """
        Create a mutant protein based on a list of mutations on a PDB file path.
        :param infile: str
        :param outfile: str the file to save the mod as.
        :param chains: lsit of str chain id in the pdb loaded.
        :return:

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol:
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            self.pymol.cmd.load(infile)
            self._chain_removal(outfile, chains)
        return 1
