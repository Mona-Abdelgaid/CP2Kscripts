import setuptools


with open('README.md', 'r') as readme:
    long_description = readme.read()

setuptools.setup(name='CP2Kscripts',
                 version="0.1.5",
                 author='Mona Abdelgaid',
                 author_email="moa59@pitt.edu",
                 description="Useful scripts for cp2k package",
                 long_description=long_description,
                 long_description_content_type='text/markdown',
                 packages=setuptools.find_packages(),
                 entry_points={'console_scripts': [
                  'interp = CP2Kscripts.interp:main',
                  'shift_coord = CP2Kscripts.shift_coord:main',
                  'vib_list = CP2Kscripts.vib_list:main',
                  'thermo = CP2Kscripts.thermo:main',
                  'gibbs_surface_all_R = CP2Kscripts.gibbs_surface_all_R:main',
                  'gibbs_surface = CP2Kscripts.gibbs_surface:main',
                  'gibbs_molecules = CP2Kscripts.gibbs_molecules:main',
                  'ldos = CP2Kscripts.ldos:main'
                  ]},
                 url="https://github.com/Mona-Abdelgaid/CP2Kscripts",
                 python_requires='>=3.6',
                 classifiers=[
                    "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: MIT License",
                    "Operating System :: OS Independent"],
                 install_requires=['matplotlib',
                                   'numpy>=1.17.2',
                                   'ase>=3.17.0'])