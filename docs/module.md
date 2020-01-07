# michelanglo_transpiler package

## Module contents

NB. Written for python 3, not tested under 2.


#### class michelanglo_transpiler.ColorItem(value: Sequence)
Bases: `object`


### \__init__(value: Sequence)
`value` is a tuple outputed from Pymol by (n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()
such as `('bismuth', 5358, (0.6196078658103943, 0.30980393290519714, 0.7098039388656616))`


* **Parameters**

    **value** (*sequence of three** (**str**, **int**, **(**float**, **float**, **float**)**)*) – sequence of (name, PyMol index, (R, G, B)) where R, G, B is under 1.



* **Variables**

    
    * **name** (*str*) – name of color


    * **index** (*int*) – PyMOL index of color


    * **rgb** (*Sequence*) – R, G, B


    * **hex** (*str*) – hex string form of color



#### class michelanglo_transpiler.ColorSwatch(colors)
Bases: `object`


### \__getitem__(index: int)

* **Parameters**

    **index** – a pymol color index



### \__init__(colors)
ColorSwatch()._swatch is a dictionary with indicing being the pymol color number. The values are ColorItem instances.
Preloading the colors is faster than querying pymol.
`print [(n, i, cmd.get_color_tuple(i)) for n,i in cmd.get_color_indices()]` in Pymol generates a good amount, but it is not the full amount.


* **Parameters**

    **colors** – a list like [(‘white’, 0, (1.0, 1.0, 1.0))]



#### class michelanglo_transpiler.PyMolTranspiler(file=None, verbose=False, validation=False, view=None, representation=None, pdb='', skip_disabled=True, job='task', run_analysis=True, \*\*settings)
Bases: `object`

The class initialises as a blank object with settings unless the file (filename of PSE file) or view and/or reps is passed.
For views see .convert_view(view_string), which processes the output of PyMOL command set_view
For representation see .convert_reps(reps_string), which process the output of PyMOL command iterate 1UBQ, print resi, resn,name,ID,reps


* **Variables**

    **swatch** (*ColorSwatch*) – all the pymol colors



### \__init__(file=None, verbose=False, validation=False, view=None, representation=None, pdb='', skip_disabled=True, job='task', run_analysis=True, \*\*settings)
Converter. `__init__` does not interact with PyMOL, so does not use the lock. Unless `run_analysis` is specified then `_postinit()` is called which does.


* **Param**

    job: this is needed for the async querying of progress in the app, but not the transpiler code itself. see .log method



* **Param**

    file: filename of PSE file.



* **Parameters**

    
    * **verbose** – print?


    * **validation** – print validation_text set for pymol?


    * **view** – the text from PymOL get_view


    * **representation** – the text from PyMOL iterate


    * **pdb** – the PDB name or code



### classmethod chain_removal_code(code, outfile, chains)
Create a mutant protein based on a list of mutations on a PDB code.


* **Parameters**

    
    * **code** – str pdb code.


    * **outfile** – str the file to save the mod as.


    * **chains** – list of str chain id in the pdb loaded.



* **Returns**

    


### classmethod chain_removal_file(infile, outfile, chains)
Create a mutant protein based on a list of mutations on a PDB file path.
:param infile: str
:param outfile: str the file to save the mod as.
:param chains: lsit of str chain id in the pdb loaded.
:return:


### static collapse_list(l: Sequence)
Given a list of residues makes a list of hyphen range string


### convert_color(uniform_non_carbon=False, inner_tabbed=1, \*\*settings)
determine what colors we have.
`{'carbon':carboncolorset,'non-carbon': colorset}`


### classmethod convert_mesh(fh, scale=0, centroid_mode='unaltered', origin=None)
Given a fh or iterable of strings, return a mesh, with optional transformations.
Note color will be lost.
Only accepts trianglular meshes!


* **Parameters**

    
    * **fh** – file handle


    * **scale** – 0 do nothing. else Angstrom size


    * **centroid_mode** – unaltered | origin | center


    * **origin** – if centroid_mode is origin get given a 3d vector.



* **Returns**

    {‘o_name’: object_name, ‘triangles’: mesh triangles}



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



### current_task( = '[2019-12-10 15:37:01.119582] idle')

### describe()
determine how and what the chains are labelled and what are their ranges.
`{'peptide': [f'{first_resi}-{last_resi}:{chain}', ..], 'hetero': [f'[{resn}]{resi}:{chain}', ..]}`


* **Return type**

    dict



### fix_structure()
Fix any issues with structure. see pymol_model_chain_segi.md for more.
empty chain issue.


### classmethod get_atom_id_of_coords(coord)
Returns the pymol atom object correspondng to coord. “Needed” for distance.


* **Parameters**

    **coord** – [x, y, z] vector



* **Returns**

    atom



### get_html(ngl='https://cdn.rawgit.com/arose/ngl/v0.10.4-1/dist/ngl.js', \*\*settings)
Returns a string to be copy-pasted into HTML code.


* **Parameters**

    
    * **ngl** – (optional) the address to ngl.js. If unspecified it gets it from the RawGit CDN


    * **viewport** – (optional) the id of the viewport div, without the hash.


    * **image** – (optional) advanced mode with clickable image?



* **Returns**

    a string.



### get_js(\*\*settings)

### get_loadfun_js(\*\*settings)

### get_reps(inner_tabbed=1, stick='sym_licorice', \*\*settings)
This method is not used.


### get_view(output='matrix', \*\*settings)
If the output is set to string, the string will be a JS command that will require the object stage to exist.


* **Parameters**

    **output** – ‘matrix’ | ‘string’



* **Returns**

    np 4x4 matrix or a NGL string



### classmethod load_pdb(file)
Loads a pdb file into a transpiler obj. and fixes it.


* **Parameters**

    **file** – str file name



* **Returns**

    self



### classmethod log(msg)

### classmethod mutate_code(code, outfile, mutations, chain)
Create a mutant protein based on a list of mutations on a PDB code.


* **Parameters**

    
    * **code** – str pdb code.


    * **outfile** – str the file to save the mod as.


    * **mutations** – list of string in the single letter format (A234P) without “p.”.


    * **chain** – str chain id in the pdb loaded.



* **Returns**

    


### classmethod mutate_file(infile, outfile, mutations, chain)
Create a mutant protein based on a list of mutations on a PDB file path.


* **Parameters**

    
    * **infile** – str


    * **outfile** – str the file to save the mod as.


    * **mutations** – list of string in the single letter format (A234P) without “p.”.


    * **chain** – str chain id in the pdb loaded.



* **Returns**

    


### parse_ss(data=None, \*\*settings)
Secondary structure


### classmethod renumber(pdb, definitions)
Fetches a pdb file into a transpiler obj.


* **Parameters**

    
    * **file** – str file name


    * **definitions** – Structure.chain_definitions


[{‘chain’: ‘A’, ‘uniprot’: ‘Q9BZ29’, ‘x’: 1605, ‘y’: 2069, ‘offset’: 1604, ‘range’: ‘1605-2069’, ‘name’: None, ‘description’: None},
:return: self


### classmethod sdf_to_pdb(infile: str, reffile: str)
A special class method to convert a sdf to pdb but with the atom index shifted so that the pdb can be cat’ed.


* **Parameters**

    
    * **infile** – sdf file


    * **reffile** – pdb file for the indices.



* **Returns**

    PDB block



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


* **Variables**

    **lock** (*threading.Lock*) – the lock. A class attribute.



### \__init__(fun)
Initialize self.  See help(type(self)) for accurate signature.


### clean_up()
Reset the Pymol instance without calling reintialise


### close_up()
Calls `clean_up` and releases the lock.


### lock( = <unlocked _thread.lock object>)

### start_up()
Starts the task in `self.fun` and takes the lock or waits.


#### michelanglo_transpiler.file_test()

#### michelanglo_transpiler.new_template_testing()

#### michelanglo_transpiler.test()