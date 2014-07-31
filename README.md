MBPol plugin for OpenMM
=======================

`mbpol` is a plugin for the `OpenMM` toolkit for molecular simulations.

It implements the `MBPol` force field, available components are:

* `MBPolReferenceElectrostaticsForce`
* `MBPolReferenceOneBodyForce`
* `MBPolReferenceTwoBodyForce`
* `MBPolReferenceThreeBodyForce`
* `MBPolReferenceDispersionForce`

The parameters of each component are defined in [`python/mbpol.xml`](https://github.com/sdsc/mbpol_openmm_plugin/blob/master/python/mbpol.xml).
As of version `0.4.0`, only the `Reference` platform, i.e. single threaded C++ on CPU, is supported.

## Requirements

* `OpenMM` built from a July 2014 or later snapshot of the Github repository <https://github.com/SimTk/openmm> master branch

## How to install

* Clone the `mbpol` plugin source from Github:
  <https://github.com/sdsc/mbpol_openmm_plugin>
* Create the `build_mbpol` folder outside of the source folder
* Configure the build by running `ccmake -i ../mbpol_openmm_plugin` inside the `build_mbpol` folder
* Make sure to set `OPENMM_BUILD_MBPOL_PLUGIN` to `ON` and `OPENMM_BUILD_MBPOL_CUDA_LIB` to `OFF` (not supported yet), `CMAKE_INSTALL_PREFIX` should contain the path to the installed `OpenMM`, by default `/usr/local/openmm`.
* Run `make` to compile the C++ library
* Run `sudo make install` to install the `mbpol` dynamic library to the
  `openmm` folder
* Run `make PythonInstall` to install the Python wrapper, it requires
  Python and `swig`, the best is to use Anaconda
* Add the OpenMM lib folder to the dynamic libraries path, generally : `export LD_LIBRARY_PATH=/usr/local/openmm/lib:$LD_LIBRARY_PATH`

## Example simulation

Simulation of a cluster of 14 water molecules:

* [IPython Notebook example on `nbviewer`](http://nbviewer.ipython.org/gist/zonca/54c7040c1cf3f583930f)
* IPython Notebook example source: <https://github.com/sdsc/mbpol_openmm_plugin/blob/master/python/water14.ipynb>
* Python Script: <https://github.com/sdsc/mbpol_openmm_plugin/blob/master/python/water14.py>
* C++ example:
  <https://github.com/sdsc/mbpol_openmm_plugin/blob/master/examples/Simulate14WaterCluster/Simulate14WaterCluster.cpp>, see the documentation in the `examples/` folder.
