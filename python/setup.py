from distutils.core import setup
from distutils.extension import Extension
import os
import sys

openmm_dir = '@OPENMM_DIR@'
mbpolplugin_header_dir = '@MBPOLPLUGIN_HEADER_DIR@'
mbpolplugin_library_dir = '@MBPOLPLUGIN_LIBRARY_DIR@'

extension = Extension(name='_mbpolplugin',
                      sources=['MBPolPluginWrapper.cpp'],
                      libraries=['OpenMM', 'MBPolPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), mbpolplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), mbpolplugin_library_dir]
                     )

setup(name='mbpolplugin',
      version='1.0',
      py_modules=['mbpolplugin'],
      ext_modules=[extension],
     )
