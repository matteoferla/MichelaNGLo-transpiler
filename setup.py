import setuptools
##from setuptools.command.install import install

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="michelanglo_transpiler",
    version="1.0.0",
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
    ],install_requires=open('requirements.txt')
)

'''
class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        ###
        install.run(self)

'''
