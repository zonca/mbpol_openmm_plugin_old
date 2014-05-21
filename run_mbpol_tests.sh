#!/usr/bin/env bash

# Copy this script to the build folder

make -j 8 
./TestReferenceAmoebaMultipoleForce
./TestReferenceAmoebaStretchBendForce
./TestReferenceMBPolTwoBodyForce
./TestReferenceMBPolThreeBodyForce
./TestReferenceMBPolDispersionForce
./TestReferenceMBPolIntegrationTest
