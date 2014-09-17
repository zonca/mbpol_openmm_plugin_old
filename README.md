MBPol plugin for OpenMM
=======================

`mbpol` is a plugin for the `OpenMM` toolkit for molecular simulations.

It implements the `MBPol` force field, available components are:

* `MBPolReferenceElectrostaticsForce`
* `MBPolReferenceOneBodyForce`
* `MBPolReferenceTwoBodyForce`
* `MBPolReferenceThreeBodyForce`
* `MBPolReferenceDispersionForce`

The parameters of each component are defined in [`python/mbpol.xml`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/mbpol.xml).
As of version `0.4.0`, only the `Reference` platform, i.e. single threaded C++ on CPU, is supported.

## How to install

See the [`INSTALL.md`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/INSTALL.md) file in the package.

## Example simulation

Simulation of a cluster of 14 water molecules:

* [IPython Notebook example on `nbviewer`](http://nbviewer.ipython.org/gist/zonca/54c7040c1cf3f583930f)
* IPython Notebook example source: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.ipynb>
* Python Script: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.py>
* C++ example:
  <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/examples/Simulate14WaterCluster/Simulate14WaterCluster.cpp>, see the documentation in the `examples/` folder.
