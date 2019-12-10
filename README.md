# MichelaNGLo-transpiler
The transpiler part of Michelaɴɢʟo which converts PyMOL to NGL

# michelanglo_transpiler package

## Module contents

NB. Written for python 3, not tested under 2. See readme.md


#### class michelanglo_transpiler.ColorItem(value: Sequence)
Bases: `object`


### \__init__(value: Sequence)
`value` is a tuple outputed from Pymol by (n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()
such as `('bismuth', 5358, (0.6196078658103943, 0.30980393290519714, 0.7098039388656616))`
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


#### class michelanglo_transpiler.ColorSwatch(colors)
Bases: `object`


### \__getitem__(index: int)

* **Parameters**

    **index** – a pymol color index



### \__init__(colors)
ColorSwatch()._swatch is a dictionary with indicing being the pymol color number. The values are ColorItem instances.
Preloading the colors is faster than querying pymol.
`print [(n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()]` in Pymol generates a good amount, but it is not the full amount.
:param colors: a list like [(‘white’, 0, (1.0, 1.0, 1.0))]


#### class michelanglo_transpiler.PyMolTranspiler(file=None, verbose=False, validation=False, view=None, representation=None, pdb='', skip_disabled=True, job='task', run_analysis=True, \*\*settings)
Bases: `object`

The class initialises as a blank object with settings unless the file (filename of PSE file) or view and/or reps is passed.
For views see .convert_view(view_string), which processes the output of PyMOL command set_view
For representation see .convert_reps(reps_string), which process the output of PyMOL command iterate 1UBQ, print resi, resn,name,ID,reps
:var swatch: all the pymol colors
:vartype swatch: ColorSwatch


### \__init__(file=None, verbose=False, validation=False, view=None, representation=None, pdb='', skip_disabled=True, job='task', run_analysis=True, \*\*settings)
Converter. `__init__` does not interact with PyMOL, so does not use the lock. Unless `run_analysis` is specified then `_postinit()` is called which does.
:param: job: this is needed for the async querying of progress in the app, but not the transpiler code itself. see .log method
:param: file: filename of PSE file.
:param verbose: print?
:param validation: print validation_text set for pymol?
:param view: the text from PymOL get_view
:param representation: the text from PyMOL iterate
:param pdb: the PDB name or code


### classmethod chain_removal_code()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### classmethod chain_removal_file()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### static collapse_list(l: Sequence)
Given a list of residues makes a list of hyphen range string


### convert_color(uniform_non_carbon=False, inner_tabbed=1, \*\*settings)
determine what colors we have.
`{'carbon':carboncolorset,'non-carbon': colorset}`
self.elemental_mapping = {}


### classmethod convert_mesh(fh, scale=0, centroid_mode='unaltered', origin=None)
Given a fh or iterable of strings, return a mesh, with optional transformations.
Note color will be lost.
Only accepts trianglular meshes!
:param fh: file handle
:param scale: 0 do nothing. else Angstrom size
:param centroid_mode: unaltered | origin | center
:param origin: if centroid_mode is origin get given a 3d vector.
:return: {‘o_name’: object_name, ‘triangles’: mesh triangles}


### convert_representation(represenation, \*\*settings)
iterate all, ID, segi, chain,resi, resn,name, elem,reps, color, ss
reps seems to be a binary number. controlling the following
\* 0th bit: sticks
\* 7th bit: line
\* 5th bit: cartoon
\* 2th bit: surface


### convert_view(view, \*\*settings)
Converts a Pymol get_view output to a NGL M4 matrix.
If the output is set to string, the string will be a JS command that will require the object stage to exist.

fog and alpha not implemented.
pymol.cmd.get(“field_of_view”))
pymol.cmd.get(“fog_start”)


* **Parameters**

    **view** – str or tuple



* **Returns**

    np 4x4 matrix or a NGL string



### current_task( = '[2019-12-10 15:10:54.538299] idle')

### describe()
determine how and what the chains are labelled and what are their ranges.
`{'peptide': [f'{first_resi}-{last_resi}:{chain}', ..], 'hetero': [f'[{resn}]{resi}:{chain}', ..]}`


### fix_structure()
Fix any issues with structure. see pymol_model_chain_segi.md for more.
empty chain issue.
:return:


### classmethod get_atom_id_of_coords(coord)
Returns the pymol atom object correspondng to coord. “Needed” for distance.
:param coord: [x, y, z] vector
:return: atom


### get_html(ngl='https://cdn.rawgit.com/arose/ngl/v0.10.4-1/dist/ngl.js', \*\*settings)
Returns a string to be copy-pasted into HTML code.
:param ngl: (optional) the address to ngl.js. If unspecified it gets it from the RawGit CDN
:param viewport: (optional) the id of the viewport div, without the hash.
:param image: (optional) advanced mode with clickable image?
:return: a string.


### get_js(\*\*settings)

### get_loadfun_js(\*\*settings)

### get_reps(inner_tabbed=1, stick='sym_licorice', \*\*settings)
This method is not used.


### get_view(output='matrix', \*\*settings)
If the output is set to string, the string will be a JS command that will require the object stage to exist.
:param output: ‘matrix’ | ‘string’
:return: np 4x4 matrix or a NGL string


### classmethod load_pdb()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### classmethod log(msg)

### classmethod mutate_code()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### classmethod mutate_file()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### parse_ss(data=None, \*\*settings)
Secondary structure


### classmethod renumber()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### classmethod sdf_to_pdb()
Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### swatch( = <michelanglo_transpiler.ColorSwatch object>)

### tmp( = '/home/matteo/Coding/MichelaNGLo-transpiler/Sphinx-docs')

### write_hmtl(template_file='test.mako', output_file='test_generated.html', \*\*kargs)

#### class michelanglo_transpiler.PyMolTranspilerDeco(fun)
Bases: `object`

Decorator for the bound methods of PyMolTranspiler that use Pymol.
The session is shared… so only one thread at the time ought to use PyMOL.
If a session raises an error, it should be caught so everyhting is cleaned closed and the error raised for the logger.
Conor has rightfully suggested that the lock should be handled by the scheduler. I.e. a request is made and the a job is added to a queue.
Currently, each extra concurrent thread simply waits or dies if it waits too long.
:var lock: the lock. A class attribute.
:vartype lock: threading.Lock


### \__init__(fun)
Initialize self.  See help(type(self)) for accurate signature.


### clean_up()
Reset the Pymol instance without calling reintialise


### close_up()
Calls `clean_up` and releases the lock.


### lock( = <unlocked _thread.lock object>)

### start_up()
Starts the task in `self.fun` and takes the lock or waits.
