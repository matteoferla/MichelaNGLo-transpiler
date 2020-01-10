from michelanglo_transpiler import PyMolTranspiler
import os, argparse
from warnings import warn
import unittest

warn('These test were used during development of the parts and teh code has changed since.')

def test_pymol_output():
    raise DeprecationWarning
    transpiler = PyMolTranspiler(verbose=True, validation=False)
    transpiler.pdb = '1UBQ'
    view = ''
    reps = ''
    data = open(os.path.join('michelanglo_app','static','pymol_demo.txt')).read().split('PyMOL>')
    for block in data:
        if 'get_view' in block:
            view = block
        elif 'iterate' in block:  # strickly lowercase as it ends in *I*terate
            reps = block
        elif not block:
            pass #empty line.
        else:
            warn('Unknown block: '+block)
    transpiler.convert_view(view)
    transpiler.convert_representation(reps)
    code = transpiler.get_html(tabbed=0)  # ngl='ngl.js'
    return transpiler



def test_new_template():
    #transpiler = PyMolTranspiler(file='git_docs/images/1ubq.pse')
    transpiler = PyMolTranspiler().transpile(file='test/1ubq.pse')
    transpiler.pdb = '1UBQ'
    transpiler.m4_alt = None
    code=transpiler.get_js(toggle_fx=True, viewport='viewport', variants=[], save_button='save_button', backgroundColor='white',tag_wrapped=True)
    transpiler.write_hmtl(template_file='test2.mako', output_file='example.html', code=code)

class Tests(unittest.TestCase):

    def test_transpile(self):
        """
        Tests the pse transpiler
        """
        transpiler = PyMolTranspiler().transpile(file='test/1ubq.pse')
        transpiler.raw_pdb

    def test_transpile2(self):
        """
        Tests the pse transpiler
        """
        settings =  {'viewport': 'viewport', 'image': None, 'uniform_non_carbon': True, 'verbose': False, 'validation': True,
     'stick_format': 'sym_licorice', 'save': True, 'backgroundcolor': 'white', 'location_viewport': 'left',
     'columns_viewport': 9, 'columns_text': 3}
        PyMolTranspiler().transpile(file='../app/michelanglo_app/demo/F.pse', **settings)



    def test_get_reps(self):
        """
        Testing the get view/reps
        """
        transpiler = PyMolTranspiler().transpile(file='test/1ubq.pse')
        v = transpiler.get_view()
        r = transpiler.get_reps()
        #print(v)
        self.assertTrue(len(v) > 0)
        #print(r)
        self.assertTrue(len(r) > 0)


if __name__ == "__main__":
    PyMolTranspiler.verbose = True
    unittest.main(verbosity=2)

