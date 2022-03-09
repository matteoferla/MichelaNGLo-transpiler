import setuptools
##from setuptools.command.install import install

import os
this_directory = os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(this_directory,'requirements.txt'), 'r') as fh:
        requirements = fh.read().strip().split()
    with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = fh.read()
except Exception:  # weird file
    requirements = ['biopython', 'pymol']
    long_description = ''

setuptools.setup(
    name="michelanglo_transpiler",
    version="1.0.1",
    author="Matteo Ferla",
    author_email="matteo@well.ox.ac.uk",
    description="Transpiler module for Michelanglo",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/matteoferla/MichelaNGLo-transpiler",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],install_requires=requirements
)

'''
class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        ###
        install.run(self)

'''
