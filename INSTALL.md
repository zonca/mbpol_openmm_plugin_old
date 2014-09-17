## Requirements

* `OpenMM` version 6.1 or later, either installed from the binary distribution or from source.

## How to install

* Clone the `mbpol` plugin source from Github:
  <https://github.com/paesanilab/mbpol_openmm_plugin>
* Create the `build_mbpol` folder outside of the source folder
* Configure the build by running `ccmake -i ../mbpol_openmm_plugin` inside the `build_mbpol` folder and press `c`
* Set:
  * `OPENMM_BUILD_MBPOL_CUDA_LIB`  `OFF` (not supported yet)
  * `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` should contain the path to the installed `OpenMM`, by default `/usr/local/openmm`.
  * `CMAKE_BUILD_TYPE` `Release`
* Press `g` to generate the configuration and exit
* Run `make` to compile the C++ library
* Run `(sudo) make install` to install the `mbpol` dynamic library to the
  `openmm` folder
* Run `make PythonInstall` to install the Python wrapper, it requires
  Python and `swig`, the best is to use Anaconda
* Add the OpenMM lib folder to the dynamic libraries path, generally : `export LD_LIBRARY_PATH=/usr/local/openmm/lib:$LD_LIBRARY_PATH`

