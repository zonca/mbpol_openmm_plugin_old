MBPol plugin for OpenMM
=======================

`mbpol` is a plugin for the `OpenMM` toolkt for molecular simulations.

It implements the `MBPol` force field, available components are:

* `MBPolReferenceElectrostaticsForce`
* `MBPolReferenceOneBodyForce`
* `MBPolReferenceTwoBodyForce`
* `MBPolReferenceThreeBodyForce`
* `MBPolReferenceDispersionForce`

As of version `0.3.0`, only the `Reference` platform, i.e. single threaded C++ on CPU, is supported.

## How to install

* Download `OpenMM` source distribution (tested 5.2 and 6.0.1)
* Clone the `mbpol` plugin source from Github to the `plugins/mbpol` folder
* Patch the `OpenMM` root `CMake` configuration file `CMakeLists.txt` by running the `apply_CMakeLists_mbpol_plugin_patch.sh` script
* Follow the `OpenMM` build instructions, making sure to set `OPENMM_BUILD_MBPOL_PLUGIN` to `ON` and `OPENMM_BUILD_MBPOL_CUDA_LIB` to `OFF` (not supported yet)
