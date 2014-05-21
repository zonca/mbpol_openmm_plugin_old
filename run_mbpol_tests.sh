#!/usr/bin/env bash

# Copy this script to the build folder

make -j 8 
./TestReferenceAmoebaMultipoleForce
./TestReferenceAmoebaStretchBendForce
./TestReferenceAmoebaVdwForce
./TestReferenceMBPolThreeBodyForce
./TestReferenceMBPolDispersionForce
./TestReferenceMBPolIntegrationTest
