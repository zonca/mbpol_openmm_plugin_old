/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
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
 * This tests the Reference implementation of ReferenceAmoebaMultipoleForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "AmoebaReferenceMultipoleForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
const double TOL = 1e-4;

static void setupWater3System( AmoebaMultipoleForce::NonbondedMethod nonbondedMethod,
                 AmoebaMultipoleForce::PolarizationType polarizationType,
                 double cutoff, std::string& testName, 
                 std::vector< double >& outputMultipoleMoments,
                 std::vector< Vec3 >& inputGrid,
                 std::vector< double >& outputGridPotential, FILE* log ){

    // beginning of Multipole setup

    System system;

    // box dimensions

    // double boxDimension                               = 1.8643;
    // Vec3 a( boxDimension, 0.0, 0.0 );
    // Vec3 b( 0.0, boxDimension, 0.0 );
    // Vec3 c( 0.0, 0.0, boxDimension );
    // system.setDefaultPeriodicBoxVectors( a, b, c );

    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 9;
    amoebaMultipoleForce->setNonbondedMethod( nonbondedMethod );
    //amoebaMultipoleForce->setPolarizationType( polarizationType );
    //amoebaMultipoleForce->setCutoffDistance( cutoff );
    //amoebaMultipoleForce->setMutualInducedTargetEpsilon( 1.0e-06 );
    //amoebaMultipoleForce->setMutualInducedMaxIterations( 500 );
    //amoebaMultipoleForce->setAEwald( 5.4459052e+00 );
    //amoebaMultipoleForce->setEwaldErrorTolerance( 1.0e-04 );

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
    }

    std::vector<double> oxygenMolecularDipole(3);
    std::vector<double> oxygenMolecularQuadrupole(9);

    oxygenMolecularDipole[0]     =   0.0000000e+00;
    oxygenMolecularDipole[1]     =   0.0000000e+00;
    oxygenMolecularDipole[2]     =   7.5561214e-03;

    oxygenMolecularQuadrupole[0] =   0.0000000e+00;
    oxygenMolecularQuadrupole[1] =   0.0000000e+00;
    oxygenMolecularQuadrupole[2] =   0.0000000e+00;
    oxygenMolecularQuadrupole[3] =   0.0000000e+00;
    oxygenMolecularQuadrupole[4] =   0.0000000e+00;
    oxygenMolecularQuadrupole[5] =   0.0000000e+00;
    oxygenMolecularQuadrupole[6] =   0.0000000e+00;
    oxygenMolecularQuadrupole[7] =   0.0000000e+00;
    oxygenMolecularQuadrupole[8] =   0.0000000e+00;

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -2.0420949e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -3.0787530e-03;

    hydrogenMolecularQuadrupole[0] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   0.0000000e+00;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
        amoebaMultipoleForce->addMultipole( -5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, jj+1, jj+2, -1,
                                            3.9000000e-01, 3.0698765e-01, 8.3700000e-04 );
        amoebaMultipoleForce->addMultipole(  2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+2, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04 );
        amoebaMultipoleForce->addMultipole(  2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+1, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04 );
    }

    // CovalentMaps

    //std::vector< int > covalentMap;
    //for( unsigned int jj = 0; jj < numberOfParticles; jj += 3 ){
    //    covalentMap.resize(0);
    //    covalentMap.push_back( jj+1 );
    //    covalentMap.push_back( jj+2 );
    //    amoebaMultipoleForce->setCovalentMap( jj, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );

    //    covalentMap.resize(0);
    //    covalentMap.push_back( jj );
    //    covalentMap.push_back( jj+1 );
    //    covalentMap.push_back( jj+2 );
    //    amoebaMultipoleForce->setCovalentMap( jj,   static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );
    //    amoebaMultipoleForce->setCovalentMap( jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );
    //    amoebaMultipoleForce->setCovalentMap( jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap );
    //
    //    covalentMap.resize(0);
    //    covalentMap.push_back( jj );
    //    amoebaMultipoleForce->setCovalentMap( jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );
    //    amoebaMultipoleForce->setCovalentMap( jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap );
    //
    //    covalentMap.resize(0);
    //    covalentMap.push_back( jj+2 );
    //    amoebaMultipoleForce->setCovalentMap( jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );
    //
    //    covalentMap.resize(0);
    //    covalentMap.push_back( jj+1 );
    //    amoebaMultipoleForce->setCovalentMap( jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap );
    //
    //} 
    system.addForce(amoebaMultipoleForce);
 
    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
    positions[3]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[4]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[5]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );
    positions[6]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[7]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[8]             = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );

    std::string platformName;
    platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);

    State state                      = context.getState(State::Forces);
    return;
}

static void testGetAndScaleInverseRs( FILE* log ) {

    std::string testName      = "testGetAndScaleInverseRs";

    
    // getAndScaleInverseRs is protected, so we need to create a wrapping class for testing it.
    class WrappedAmoebaReferenceMultipoleForce : public AmoebaReferenceMultipoleForce {
        public:
        int wrapGetAndScaleInverseRs(
                RealOpenMM dampI, RealOpenMM dampJ,
                RealOpenMM tholeI, RealOpenMM tholeJ,
                RealOpenMM r, bool justScale, RealOpenMM & damp, std::vector<RealOpenMM>& rrI
                )   { 
                        getAndScaleInverseRs(dampI, dampJ, tholeI, tholeJ, r, justScale, damp, rrI);
                    }
    };

    RealOpenMM damp=10.;
    RealOpenMM dampO=0.306988;
    RealOpenMM dampH=0.28135;
    std::vector<RealOpenMM> rrI(4);
    RealOpenMM r=9.860634018e-01; // from Water3 test
    RealOpenMM thole=0.3900;

    WrappedAmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = new WrappedAmoebaReferenceMultipoleForce();;
    amoebaReferenceMultipoleForce->wrapGetAndScaleInverseRs( dampO, dampH,
                          thole, thole, r, false, damp, rrI);

    //ASSERT_EQUAL_TOL_MOD(0., rrI[0], 1e-5, testName);
    //ASSERT_EQUAL_TOL_MOD(5.324612470e-01, rrI[1], 1e-5, testName);
    //ASSERT_EQUAL_TOL_MOD(4.747626558e-02, rrI[2], 1e-5, testName);
    //ASSERT_EQUAL_TOL_MOD(             0., rrI[3], 1e-5, testName);

}

static void testWater3Distances( FILE* log ) {

    std::string testName      = "testWaterDistances";
    
    int numberOfParticles     = 9;
    double cutoff             = 0.70;

    std::vector<Vec3> forces;
    double energy;
    std::vector<double> outputMultipoleMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;

    setupWater3System( AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Mutual, 
                          cutoff, testName,
                       outputMultipoleMoments, inputGrid, outputGridPotential, log );

    //State state                      = water3Context->getState(State::Forces);
    // forces                           = state.getForces();

    //std::vector<RealOpenMM> rrI(4);
    //getAndScaleInverseRs( particleI.dampingFactor, particleK.dampingFactor,
    //                      particleI.thole, particleK.thole, r, false, damp, rrI );

//    std::vector<double> tinkerMoments(4);
//
//    tinkerMoments[0]  =   0.0000000e+00;
//    tinkerMoments[1]  =  -9.1118361e+00;
//    tinkerMoments[2]  =   2.8371876e+01;
//    tinkerMoments[3]  =   5.1518898e+01;
////    tinkerMoments[4]  =  -1.0768808e-01; // Quadrupole moments are not uniquely defined when using periodic boundary conditions
////    tinkerMoments[5]  =  -9.0458124e-01;
////    tinkerMoments[6]  =   1.8460385e+00;
////    tinkerMoments[7]  =  -9.0458124e-01;
////    tinkerMoments[8]  =   6.4395591e-02;
////    tinkerMoments[9]  =   1.6692567e-01;
////    tinkerMoments[10] =   1.8460385e-00;
////    tinkerMoments[11] =   1.6692567e-01;
////    tinkerMoments[12] =   4.3292490e-02;
//
//    double tolerance = 1.0e-04;
//    if( log ){
//        (void) fprintf( log, "%s RelativeDifference Tinker OpenMM\n", testName.c_str() );
//    }
//    for( unsigned int ii = 0; ii < tinkerMoments.size(); ii++ ){
//        double difference;
//        if( fabs( tinkerMoments[ii] ) > 0.0 ){
//            difference = fabs( outputMultipoleMoments[ii] - tinkerMoments[ii] )/fabs( tinkerMoments[ii] );
//        } else {
//            difference = fabs( outputMultipoleMoments[ii] - tinkerMoments[ii] );
//        }
//        if( log ){
//            (void) fprintf( log, "%2d %15.7e %15.7e %15.7e\n", ii, difference, tinkerMoments[ii], outputMultipoleMoments[ii] );
//        }
//
//        if( difference > tolerance ){
//            std::stringstream details;
//            details << testName << "Multipole moment " << ii << " does not agree w/ TINKER computed moments: OpenMM=" << outputMultipoleMoments[ii];
//            details << " TINKER=" <<  tinkerMoments[ii]  << " difference=" << difference;
//            throwException(__FILE__, __LINE__, details.str());
//        }
//
//    }

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceAmoebaMultipoleForce running test..." << std::endl;

        FILE* log = NULL;

        testGetAndScaleInverseRs( log );

        // water 3 mbpol
        testWater3Distances( log );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
