# MichelaNGLo-transpiler
The transpiler part of Michelaɴɢʟo which converts PyMOL to NGL.

The site [michelanglo.sgc.ox.ac.uk](https://michelanglo.sgc.ox.ac.uk) depends on three repos:
* [MichelaNGLo-app](https://github.com/matteoferla/MichelaNGLo)
* **MichelaNGLo-transpiler**
* [MichelaNGLo-protein-module](https://github.com/matteoferla/MichelaNGLo-protein-module)

The top-level repo is [Michelanglo-and-Venus](https://github.com/matteoferla/Michelanglo-and-Venus)

This module is self-standing, so can be used outside of the app or protein-module.

## Implementation

The conversion file, `michelanglo/transpiler.py` and the files in `michelanglo/transpiler_templates` are all that is required to use locally.
They are in `michelanglo_app` to avoid allowing the app to do a relative import beyond the top-level package (`michelanglo_app`).

For notes about the details about the conversion see [conversion.md](docs/conversion.md)
For notes about the transpiler see [transpiler.md](docs/transpiler.md).

## Further details

* [Python package details](docs/module.md)
* [conversion details](docs/conversion.md)
* [Note on conversion](docs/notes_on_view_conversion.md)
* [PyMol models and segi](docs/PyMOL_model_chains_segi.md)
* [Note on surfaces](docs/notes_on_surfaces.md)
* [Note on mmCIF and mmTF](docs/notes_on_mmCIF.md)
* [Multiple Pymol sessions](docs/pymol_sessions.md)