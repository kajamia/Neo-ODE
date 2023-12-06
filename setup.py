from skbuild import setup

import sys
_cmake_args = []

if 'win32' in sys.platform:
    _cmake_args += ['-G', 'MinGW Makefiles']

    
setup(
    name="Neoode",
    version="0.0.2",
    author="Kaja M., Vivien R. Dániel L., Norbert H.",
    license="MIT",
    packages=["Neoode"],
    cmake_args=_cmake_args
)
