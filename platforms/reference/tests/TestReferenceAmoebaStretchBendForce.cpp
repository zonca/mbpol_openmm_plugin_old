/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the Reference implementation of ReferenceAmoebaStretchBendForce.
 */

#include "openmm/internal/AssertionUtilities.h"
//#include "AmoebaTinkerParameterFile.h"
const double DegreesToRadians = 3.14159265/180.0;
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

const double TOL = 1e-4;
#define PI_M               3.141592653589
#define RADIAN            57.29577951308

void testOneStretchBend( FILE* log ) {

    System system;
    int numberOfParticles = 3;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaStretchBendForce* amoebaStretchBendForce = new AmoebaStretchBendForce();

    double abLength         = 0.144800000E+01;
    double cbLength         = 0.101500000E+01;
    double angleStretchBend = 0.108500000E+03*DegreesToRadians;
    //double kStretchBend     = 0.750491578E-01;
    double kStretchBend     = 1.0;

    amoebaStretchBendForce->addStretchBend(0, 1, 2);

    system.addForce(amoebaStretchBendForce);
    Context context(system, integrator, Platform::getPlatformByName( "Reference"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );

    context.setPositions(positions);
    
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();
    double cal2joule = 4.184;

    std::vector<Vec3> expectedForces(numberOfParticles);
    expectedForces[0]         = Vec3( -1.029233628e-01,  1.752006876e-01, -2.394228296e-01  );
    expectedForces[1]         = Vec3(  1.238286503e-01, -9.713944883e-02,  9.278441270e-02  );
    expectedForces[2]         = Vec3( -1.992936921e-02, -8.084103617e-02,  1.660930712e-01  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            expectedForces[i][j] *= cal2joule*10;
        }
    }
    double tolerance=1e-4;
    double expectedEnergy = -15.9939592*cal2joule;

    std::cout  << std::endl << "Forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
         std::cout << forces[i] << " Kcal/mol/A " << std::endl;
    }

    std::cout  << std::endl << "Expected forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
         std::cout << expectedForces[i] << " Kcal/mol/A " << std::endl;
    }

    double energy = state.getPotentialEnergy();

    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
    }
    ASSERT_EQUAL_TOL( expectedEnergy, energy, tolerance );

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceAmoebaStretchBendForce running test..." << std::endl;

        FILE* log = NULL;
        //FILE* log = stderr;
        //FILE* log = fopen( "AmoebaStretchBendForce1.log", "w" );;
        testOneStretchBend( log );
#ifdef AMOEBA_DEBUG
        if( log && log != stderr )
            (void) fclose( log );
#endif

    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    //std::cout << "PASS - Test succeeded." << std::endl;
    std::cout << "Done" << std::endl;
    return 0;
}
