########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo  --- IO parts plus basic cleaning of PDB block
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

########################################################################################################################

from warnings import warn
from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint
import pymol2, re

###############################################################

class  PyMolTranspiler_io:

    def load_pdb(self, file, outfile=None, mod_fx=None):
        """
        Loads a pdb file into a transpiler obj. and fixes it.
        The round trip is to prevent anything malicious being sent.

        :param file: str file name
        :return: self

        **PyMOL session**: self-contained.
        """
        with pymol2.PyMOL() as self.pymol: #pymol2.PyMOL()
            self.pymol.cmd.set('fetch_path', self.temporary_folder)
            self.pymol.cmd.load(file)
            extension = file.split('.')[-1]
            headers = []
            gather_ss = True
            if extension == 'pdb':
                with open(file) as w:
                    headers = [row.replace('"','').replace("'",'').replace("\\",'') for row in w if any([k in row for k in ('LINK', 'HELIX', 'SHEET')])]
                    if any(['HELIX', 'SHEET' in headers]):
                        gather_ss = False
                if outfile is None:
                    outfile = file
            else:
                if outfile is None:
                    outfile = '.'.join(file.split('.')[:-1])+'.pdb'
            if mod_fx:
                mod_fx()
            self.raw_pdb = self.remove_anisou(self.pymol.cmd.get_pdbstr())
            ## fix the segi and multiple object problem.
            self.fix_structure()
            ## add SS
            if gather_ss:
                myspace = {'data': []}
                # myspace['data'] is the same as self.atoms, which is "kind of the same" as pymol.cmd.get_model('..').atoms
                self.pymol.cmd.iterate('all', self._iterate_cmd, space=myspace)
                self.parse_ss(myspace['data'])
                self.raw_pdb = '\n'.join(self.ss)+'\n'+ self.raw_pdb
            else:
                self.raw_pdb = '\n'.join(headers)+'\n'+ self.raw_pdb
            return self

    @property
    def pdb_block(self):
        if self.raw_pdb == '':
            warn('raw_PDB is empty')
            self.raw_pdb = self.pymol.cmd.get_pdbstr()
        if self.headers:
            return '\n'.join(self.headers) + '\n' + self.remove_anisou(self.raw_pdb.lstrip())
        else:
            return '\n'.join(self.ss) + '\n' + self.remove_anisou(self.raw_pdb.lstrip())

    def fix_structure(self):
        """
        Fix any issues with structure. see self.pymol_model_chain_segi.md for more.
        empty chain issue.

        **PyMOL session**: dependent. Requires sigleton.
        """

        # whereas a chain can be called ?, it causes problems. So these are strictly JS \w characters.
        # Only latin-1 is okay in NGL. Any character above U+00FF will be rendered as the last two bytes. (U+01FF will be U+00FF say)
        #non-ascii are not okay in PyMol
        chaingen = self.get_new_letter()
        objs = self.pymol.cmd.get_names(enabled_only=1)
        prime_chains = self.get_chains(objs[0])
        for on in objs[1:]:
            for c in self.get_chains(on):
                if not c: # missing chain ID is still causing issues.
                    new_chain = next(chaingen)
                    self.pymol.cmd.alter(f"{on} and chain ''", f'chain="{new_chain}"')
                elif c in prime_chains:
                    new_chain = next(chaingen)
                    self.pymol.cmd.alter(f"{on} and chain {c}", f'chain="{new_chain}"')
                else:
                    prime_chains.add(c)
        # selenomethionine to methionine
        if self.pymol.cmd.select('resn MSE') > 0:
            self.pymol.cmd.alter('resn MSE and element SE', 'name=" SD "')
            self.pymol.cmd.alter('resn MSE and element SE', 'element="S"')
            self.pymol.cmd.alter('resn MSE', 'resn="MET"')
            self.pymol.cmd.sort('all')
        self.pymol.cmd.alter("all", "segi=''") # not needed. NGL does not recognise segi. Currently writtten to ignore it.
        self.pymol.cmd.sort('all')
        # The delete states shortcut does not work:
        self.pymol.cmd.create('mike_combined','enabled',1) #the 1 means that only the first "state" = model is used.
        for on in self.pymol.cmd.get_names_of_type('object:molecule'):
            if on != 'mike_combined':
                self.pymol.cmd.delete(on)

    def parse_ss(self, data=None, **settings):
        """
        PDB block Secondary structure maker
        """
        def _deal_with():
            if ss_last == 'H':  # previous was the other type
                self.ss.append('{typos}  {ss_count: >3} {ss_count: >3} {resn_start} {chain} {resi_start: >4}  {resn_end} {chain} {resi_end: >4} {h_class: >2}                                  {length: >2}'.format(
                    typos='HELIX',
                    ss_count=ss_count[ss_last],
                    resn_start=resn_start,
                    resi_start=resi_start,
                    resn_end=resn_last,
                    resi_end=resi_last,
                    chain=chain,
                    h_class=1,
                    length=int(resi_last) - int(resi_start) # already polished
                ))
                ss_count[ss_last] += 1
            elif ss_last == 'S':  # previous was the other type
                self.ss.append('{typos}  {ss_count: >3} {ss_count: >2}S 1 {resn_start} {chain}{resi_start: >4}  {resn_end} {chain}{resi_end: >4}  0'.format(
                    typos='SHEET',
                    ss_count=ss_count[ss_last],
                    resn_start=resn_start,
                    resi_start=resi_start,
                    resn_end=resn_last,
                    resi_end=resi_last,
                    chain=chain,
                    h_class=0,
                    length=resi_last - resi_start
                ))
                ss_count[ss_last] += 1

        self.ss = []
        if data is None:
            myspace = {'data': []}
            self.pymol.cmd.iterate('all', self._iterate_cmd, space=myspace)
            data = myspace['data']
        ss_last = 'L'
        resi_start = '0'
        resn_start = 'XXX'
        resi_last = '0'
        resn_last = 'XXX'
        ss_count = {'H': 1, 'S': 1, 'L': 0}
        chain = 'X'
        icode_polish = lambda resi: int(re.search(r'(\d+)', resi).group(1))
        for line in data:  # ss_list:
            if line['name'] == 'CA':
                (resi_this, ss_this, resn_this, chain) = (line['resi'], line['ss'], line['resn'], line['chain'])
                if ss_last != ss_this:
                    # deal with previous first
                    _deal_with()
                    # deal with current
                    if ss_this in ('S', 'H'):  # start of a new
                        resi_start = icode_polish(resi_this)
                        resn_start = resn_this
                        ss_last = ss_this
                # move on
                resi_last = icode_polish(resi_this)
                resn_last = resn_this
                ss_last = ss_this
        _deal_with()
        return self

    @staticmethod
    def remove_anisou(block):
        return '\n'.join([r for r in block.split('\n') if 'ANISOU' not in r])


