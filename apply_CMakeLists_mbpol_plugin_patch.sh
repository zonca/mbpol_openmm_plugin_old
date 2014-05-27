#!/usr/bin/env bash

# patch the root folder OpenMM CMakeLists to activate the mbpol plugin
patch ../../CMakeLists.txt < CMakeLists_mbpol_plugin.patch
