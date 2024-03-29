/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "MBPolReferenceKernels.h"
#include "MBPolReferenceOneBodyForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include "MBPolReferenceThreeBodyForce.h"
#include "MBPolReferenceDispersionForce.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/internal/MBPolElectrostaticsForceImpl.h"
#include "openmm/internal/MBPolTwoBodyForceImpl.h"
#include "openmm/internal/MBPolThreeBodyForceImpl.h"
#include "openmm/internal/MBPolDispersionForceImpl.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/NonbondedForceImpl.h"

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(RealVec*) data->periodicBoxSize;
}

ReferenceCalcMBPolOneBodyForceKernel::ReferenceCalcMBPolOneBodyForceKernel(std::string name, const Platform& platform, const OpenMM::System& system) :
                   CalcMBPolOneBodyForceKernel(name, platform), system(system) {
}

ReferenceCalcMBPolOneBodyForceKernel::~ReferenceCalcMBPolOneBodyForceKernel() {
}

void ReferenceCalcMBPolOneBodyForceKernel::initialize(const OpenMM::System& system, const MBPolOneBodyForce& force) {

    numOneBodys = force.getNumOneBodys();
    for ( int ii = 0; ii < numOneBodys; ii++) {
        int particle1Index, particle2Index, particle3Index;
        double lengthAB, lengthCB, angle, k;
        force.getOneBodyParameters(ii, particle1Index, particle2Index, particle3Index);
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
    }
}

double ReferenceCalcMBPolOneBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    MBPolReferenceOneBodyForce MBPolReferenceOneBodyForce;
    RealOpenMM energy      = MBPolReferenceOneBodyForce.calculateForceAndEnergy( numOneBodys, posData, particle1, particle2, particle3,
                                                                                      lengthABParameters, lengthCBParameters, angleParameters, kParameters, forceData );
    return static_cast<double>(energy);
}

void ReferenceCalcMBPolOneBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolOneBodyForce& force) {
    if (numOneBodys != force.getNumOneBodys())
        throw OpenMMException("updateParametersInContext: The number of stretch-bends has changed");

    // Record the values.

    for (int i = 0; i < numOneBodys; ++i) {
        int particle1Index, particle2Index, particle3Index;
        force.getOneBodyParameters(i, particle1Index, particle2Index, particle3Index);
        if (particle1Index != particle1[i] || particle2Index != particle2[i] || particle3Index != particle3[i])
            throw OpenMMException("updateParametersInContext: The set of particles in a stretch-bend has changed");
    }
}

/* -------------------------------------------------------------------------- *
 *                             MBPolElectrostatics                                *
 * -------------------------------------------------------------------------- */

ReferenceCalcMBPolElectrostaticsForceKernel::ReferenceCalcMBPolElectrostaticsForceKernel(std::string name, const Platform& platform, const OpenMM::System& system) : 
         CalcMBPolElectrostaticsForceKernel(name, platform), system(system), numElectrostatics(0), mutualInducedMaxIterations(200), mutualInducedTargetEpsilon(1.0e-03),
                                                         usePme(false),alphaEwald(0.0), cutoffDistance(1.0) {  

}

ReferenceCalcMBPolElectrostaticsForceKernel::~ReferenceCalcMBPolElectrostaticsForceKernel() {
}

void ReferenceCalcMBPolElectrostaticsForceKernel::initialize(const OpenMM::System& system, const MBPolElectrostaticsForce& force) {

    numElectrostatics   = force.getNumElectrostatics();

    charges.resize(numElectrostatics);
    dipoles.resize(3*numElectrostatics);
    quadrupoles.resize(9*numElectrostatics);
    tholes.resize(5*numElectrostatics);
    dampingFactors.resize(numElectrostatics);
    polarity.resize(numElectrostatics);
    axisTypes.resize(numElectrostatics);
    multipoleAtomZs.resize(numElectrostatics);
    multipoleAtomXs.resize(numElectrostatics);
    multipoleAtomYs.resize(numElectrostatics);
    multipoleAtomCovalentInfo.resize(numElectrostatics);

    int dipoleIndex      = 0;
    int quadrupoleIndex  = 0;
    int tholeIndex       = 0;
    int maxCovalentRange = 0;
    double totalCharge   = 0.0;
    for( int ii = 0; ii < numElectrostatics; ii++ ){

        // multipoles

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
        std::vector<double> quadrupolesD;
        std::vector<double> tholesD;
        force.getElectrostaticsParameters(ii, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                     tholesD, dampingFactorD, polarityD );

        totalCharge                       += charge;
        axisTypes[ii]                      = axisType;
        multipoleAtomZs[ii]                = multipoleAtomZ;
        multipoleAtomXs[ii]                = multipoleAtomX;
        multipoleAtomYs[ii]                = multipoleAtomY;

        charges[ii]                        = static_cast<RealOpenMM>(charge);

        dampingFactors[ii]                 = static_cast<RealOpenMM>(dampingFactorD);
        polarity[ii]                       = static_cast<RealOpenMM>(polarityD);

        dipoles[dipoleIndex++]             = static_cast<RealOpenMM>(dipolesD[0]);
        dipoles[dipoleIndex++]             = static_cast<RealOpenMM>(dipolesD[1]);
        dipoles[dipoleIndex++]             = static_cast<RealOpenMM>(dipolesD[2]);
        
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[0]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[1]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[2]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[3]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[4]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[5]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[6]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[7]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[8]);

        tholes[tholeIndex++]             = static_cast<RealOpenMM>(tholesD[0]);
        tholes[tholeIndex++]             = static_cast<RealOpenMM>(tholesD[1]);
        tholes[tholeIndex++]             = static_cast<RealOpenMM>(tholesD[2]);
        tholes[tholeIndex++]             = static_cast<RealOpenMM>(tholesD[3]);
        tholes[tholeIndex++]             = static_cast<RealOpenMM>(tholesD[4]);

        // covalent info

        std::vector< std::vector<int> > covalentLists;
        force.getCovalentMaps(ii, covalentLists );
        multipoleAtomCovalentInfo[ii] = covalentLists;

    }

    polarizationType = force.getPolarizationType();
    if( polarizationType == MBPolElectrostaticsForce::Mutual ){
        mutualInducedMaxIterations = force.getMutualInducedMaxIterations();
        mutualInducedTargetEpsilon = force.getMutualInducedTargetEpsilon();
    }

    includeChargeRedistribution = force.getIncludeChargeRedistribution();

    // PME

    nonbondedMethod  = force.getNonbondedMethod();
    if( nonbondedMethod == MBPolElectrostaticsForce::PME ){
        usePme     = true;
        alphaEwald = force.getAEwald();
        cutoffDistance = force.getCutoffDistance();
        force.getPmeGridDimensions(pmeGridDimension);
        if (pmeGridDimension[0] == 0 || alphaEwald == 0.0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            int gridSizeX, gridSizeY, gridSizeZ;
            NonbondedForceImpl::calcPMEParameters(system, nb, alphaEwald, gridSizeX, gridSizeY, gridSizeZ);
            pmeGridDimension[0] = gridSizeX;
            pmeGridDimension[1] = gridSizeY;
            pmeGridDimension[2] = gridSizeZ;
        }    
    } else {
        usePme = false;
    }
    return;
}

MBPolReferenceElectrostaticsForce* ReferenceCalcMBPolElectrostaticsForceKernel::setupMBPolReferenceElectrostaticsForce(ContextImpl& context )
{

    // mbpolReferenceElectrostaticsForce is set to MBPolReferenceGeneralizedKirkwoodForce if MBPolGeneralizedKirkwoodForce is present
    // mbpolReferenceElectrostaticsForce is set to MBPolReferencePmeElectrostaticsForce if 'usePme' is set
    // mbpolReferenceElectrostaticsForce is set to MBPolReferenceElectrostaticsForce otherwise

    MBPolReferenceElectrostaticsForce* mbpolReferenceElectrostaticsForce = NULL;
    if( usePme ) {

         MBPolReferencePmeElectrostaticsForce* mbpolReferencePmeElectrostaticsForce = new MBPolReferencePmeElectrostaticsForce( );
         mbpolReferencePmeElectrostaticsForce->setAlphaEwald( alphaEwald );
         mbpolReferencePmeElectrostaticsForce->setCutoffDistance( cutoffDistance );
         mbpolReferencePmeElectrostaticsForce->setPmeGridDimensions( pmeGridDimension );
         RealVec& box = extractBoxSize(context);
         double minAllowedSize = 1.999999*cutoffDistance;
         if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
         }
         mbpolReferencePmeElectrostaticsForce->setPeriodicBoxSize(box);
         mbpolReferenceElectrostaticsForce = static_cast<MBPolReferenceElectrostaticsForce*>(mbpolReferencePmeElectrostaticsForce);

    } else {
         mbpolReferenceElectrostaticsForce = new MBPolReferenceElectrostaticsForce( MBPolReferenceElectrostaticsForce::NoCutoff );
    }

    // set polarization type

    if( polarizationType == MBPolElectrostaticsForce::Mutual ){
        mbpolReferenceElectrostaticsForce->setPolarizationType( MBPolReferenceElectrostaticsForce::Mutual );
        mbpolReferenceElectrostaticsForce->setMutualInducedDipoleTargetEpsilon( mutualInducedTargetEpsilon );
        mbpolReferenceElectrostaticsForce->setMaximumMutualInducedDipoleIterations( mutualInducedMaxIterations );
    } else if( polarizationType == MBPolElectrostaticsForce::Direct ){
        mbpolReferenceElectrostaticsForce->setPolarizationType( MBPolReferenceElectrostaticsForce::Direct );
    } else {
        throw OpenMMException("Polarization type not recognzied." );
    }

    mbpolReferenceElectrostaticsForce->setIncludeChargeRedistribution(includeChargeRedistribution);

    return mbpolReferenceElectrostaticsForce;

}

double ReferenceCalcMBPolElectrostaticsForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    MBPolReferenceElectrostaticsForce* mbpolReferenceElectrostaticsForce = setupMBPolReferenceElectrostaticsForce( context );

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy          = mbpolReferenceElectrostaticsForce->calculateForceAndEnergy( posData, charges, dipoles, quadrupoles, tholes,
                                                                                         dampingFactors, polarity, axisTypes, 
                                                                                         multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                                         multipoleAtomCovalentInfo, forceData);

    delete mbpolReferenceElectrostaticsForce;

    return static_cast<double>(energy);
}

void ReferenceCalcMBPolElectrostaticsForceKernel::getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                                        std::vector< double >& outputElectrostaticPotential ){

    MBPolReferenceElectrostaticsForce* mbpolReferenceElectrostaticsForce = setupMBPolReferenceElectrostaticsForce( context );
    vector<RealVec>& posData                                     = extractPositions(context);
    vector<RealVec> grid( inputGrid.size() );
    vector<RealOpenMM> potential( inputGrid.size() );
    for( unsigned int ii = 0; ii < inputGrid.size(); ii++ ){
        grid[ii] = inputGrid[ii];
    }
    mbpolReferenceElectrostaticsForce->calculateElectrostaticPotential( posData, charges, dipoles, quadrupoles, tholes,
                                                                    dampingFactors, polarity, axisTypes, 
                                                                    multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                    multipoleAtomCovalentInfo, grid, potential );

    outputElectrostaticPotential.resize( inputGrid.size() );
    for( unsigned int ii = 0; ii < inputGrid.size(); ii++ ){
        outputElectrostaticPotential[ii] = potential[ii];
    }

    delete mbpolReferenceElectrostaticsForce;

    return;
}

void ReferenceCalcMBPolElectrostaticsForceKernel::getSystemElectrostaticsMoments(ContextImpl& context, std::vector< double >& outputElectrostaticsMoments){

    // retrieve masses

    const OpenMM::System& system             = context.getSystem();
    vector<RealOpenMM> masses;
    for (int i = 0; i <  system.getNumParticles(); ++i) {
        masses.push_back( static_cast<RealOpenMM>(system.getParticleMass(i)) );
    }    

    MBPolReferenceElectrostaticsForce* mbpolReferenceElectrostaticsForce = setupMBPolReferenceElectrostaticsForce( context );
    vector<RealVec>& posData                                     = extractPositions(context);
    mbpolReferenceElectrostaticsForce->calculateMBPolSystemElectrostaticsMoments( masses, posData, charges, dipoles, quadrupoles, tholes,
                                                                          dampingFactors, polarity, axisTypes, 
                                                                          multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                          multipoleAtomCovalentInfo, outputElectrostaticsMoments );

    delete mbpolReferenceElectrostaticsForce;

    return;
}

void ReferenceCalcMBPolElectrostaticsForceKernel::copyParametersToContext(ContextImpl& context, const MBPolElectrostaticsForce& force) {
    if (numElectrostatics != force.getNumElectrostatics())
        throw OpenMMException("updateParametersInContext: The number of multipoles has changed");

    // Record the values.

    int dipoleIndex = 0;
    int quadrupoleIndex = 0;
    int tholeIndex = 0;
    for (int i = 0; i < numElectrostatics; ++i) {
        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
        std::vector<double> tholeD;
        std::vector<double> quadrupolesD;
        force.getElectrostaticsParameters(i, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, tholeD, dampingFactorD, polarityD);
        axisTypes[i] = axisType;
        multipoleAtomZs[i] = multipoleAtomZ;
        multipoleAtomXs[i] = multipoleAtomX;
        multipoleAtomYs[i] = multipoleAtomY;
        charges[i] = (RealOpenMM) charge;
        tholes[tholeIndex++] = (RealOpenMM) tholeD[0];
        tholes[tholeIndex++] = (RealOpenMM) tholeD[1];
        tholes[tholeIndex++] = (RealOpenMM) tholeD[2];
        tholes[tholeIndex++] = (RealOpenMM) tholeD[3];
        tholes[tholeIndex++] = (RealOpenMM) tholeD[4];
        dampingFactors[i] = (RealOpenMM) dampingFactorD;
        polarity[i] = (RealOpenMM) polarityD;
        dipoles[dipoleIndex++] = (RealOpenMM) dipolesD[0];
        dipoles[dipoleIndex++] = (RealOpenMM) dipolesD[1];
        dipoles[dipoleIndex++] = (RealOpenMM) dipolesD[2];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[0];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[1];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[2];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[3];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[4];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[5];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[6];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[7];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[8];
    }
}


ReferenceCalcMBPolTwoBodyForceKernel::ReferenceCalcMBPolTwoBodyForceKernel(std::string name, const Platform& platform, const OpenMM::System& system) :
       CalcMBPolTwoBodyForceKernel(name, platform), system(system) {
    useCutoff = 0;
    usePBC = 0;
    cutoff = 1.0e+10;
    neighborList = NULL;
}

ReferenceCalcMBPolTwoBodyForceKernel::~ReferenceCalcMBPolTwoBodyForceKernel() {
    if( neighborList ){
        delete neighborList;
    } 
}

void ReferenceCalcMBPolTwoBodyForceKernel::initialize(const OpenMM::System& system, const MBPolTwoBodyForce& force) {

    // per-particle parameters

    numParticles = force.getNumMolecules();
    allParticleIndices.resize(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){

        std::vector<int> particleIndices;
        force.getParticleParameters(ii, particleIndices );
        allParticleIndices[ii] = particleIndices;

    }

    useCutoff              = (force.getNonbondedMethod() != MBPolTwoBodyForce::NoCutoff);
    usePBC                 = (force.getNonbondedMethod() == MBPolTwoBodyForce::CutoffPeriodic);
    cutoff                 = force.getCutoff();
    neighborList           = useCutoff ? new NeighborList() : NULL;

}

double ReferenceCalcMBPolTwoBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& allPosData   = extractPositions(context);
    vector<RealVec> posData;
    posData.resize(numParticles);
    vector<set<int> > allExclusions;
    allExclusions.resize(numParticles);
    // posData has only oxygens
    for( int ii = 0; ii < numParticles; ii++ ){
        posData[ii] = allPosData[allParticleIndices[ii][0]];
    }
    vector<RealVec>& forceData = extractForces(context);
    MBPolReferenceTwoBodyForce TwoBodyForce;
    RealOpenMM energy;
    TwoBodyForce.setCutoff( cutoff );
    // neighborList created only with oxygens, then allParticleIndices is used to get reference to the hydrogens
    computeNeighborListVoxelHash( *neighborList, numParticles, posData, allExclusions, extractBoxSize(context), usePBC, cutoff, 0.0, false);
    if( usePBC ){
        TwoBodyForce.setNonbondedMethod( MBPolReferenceTwoBodyForce::CutoffPeriodic);
        RealVec& box = extractBoxSize(context);
        double minAllowedSize = 1.999999*cutoff;
        if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
            throw OpenMMException("The periodic box size has decreased to less than twice the cutoff.");
        }
        TwoBodyForce.setPeriodicBox(box);
    } else {
        TwoBodyForce.setNonbondedMethod( MBPolReferenceTwoBodyForce::CutoffNonPeriodic);
    }
    // here we need allPosData, every atom!
    energy  = TwoBodyForce.calculateForceAndEnergy( numParticles, allPosData, allParticleIndices, *neighborList, forceData);

    return static_cast<double>(energy);
}

void ReferenceCalcMBPolTwoBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolTwoBodyForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {

        std::vector<int> particleIndices;
        force.getParticleParameters(i, particleIndices);
        allParticleIndices[i] = particleIndices;

    }
}

ReferenceCalcMBPolThreeBodyForceKernel::ReferenceCalcMBPolThreeBodyForceKernel(std::string name, const Platform& platform, const OpenMM::System& system) :
       CalcMBPolThreeBodyForceKernel(name, platform), system(system) {
    useCutoff = 0;
    usePBC = 0;
    cutoff = 1.0e+10;
    neighborList = NULL;
}

ReferenceCalcMBPolThreeBodyForceKernel::~ReferenceCalcMBPolThreeBodyForceKernel() {
    if( neighborList ){
        delete neighborList;
    }
}

void ReferenceCalcMBPolThreeBodyForceKernel::initialize(const OpenMM::System& system, const MBPolThreeBodyForce& force) {

    // per-particle parameters

    numParticles = force.getNumMolecules();
    allParticleIndices.resize(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){

        std::vector<int> particleIndices;
        force.getParticleParameters(ii, particleIndices );
        allParticleIndices[ii] = particleIndices;

    }

    useCutoff              = (force.getNonbondedMethod() != MBPolThreeBodyForce::NoCutoff);
    usePBC                 = (force.getNonbondedMethod() == MBPolThreeBodyForce::CutoffPeriodic);
    cutoff                 = force.getCutoff();
    neighborList           = useCutoff ? new ThreeNeighborList() : NULL;

}

double ReferenceCalcMBPolThreeBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& allPosData   = extractPositions(context);
    vector<RealVec> posData;
    posData.resize(numParticles);
    vector<set<int> > allExclusions;
    allExclusions.resize(numParticles);
    // posData has only oxygens
    for( int ii = 0; ii < numParticles; ii++ ){
        posData[ii] = allPosData[allParticleIndices[ii][0]];
    }
    vector<RealVec>& forceData = extractForces(context);
    MBPolReferenceThreeBodyForce force;
    RealOpenMM energy;
    force.setCutoff( cutoff );
    // neighborList created only with oxygens, then allParticleIndices is used to get reference to the hydrogens
    computeThreeNeighborListVoxelHash( *neighborList, numParticles, posData, extractBoxSize(context), usePBC, cutoff, 0.0);
    if( usePBC ){
        force.setNonbondedMethod( MBPolReferenceThreeBodyForce::CutoffPeriodic);
        RealVec& box = extractBoxSize(context);
        double minAllowedSize = 1.999999*cutoff;
        if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
            throw OpenMMException("The periodic box size has decreased to less than twice the cutoff.");
        }
        force.setPeriodicBox(box);
    } else {
        force.setNonbondedMethod( MBPolReferenceThreeBodyForce::CutoffNonPeriodic);
    }
    // here we need allPosData, every atom!
    energy  = force.calculateForceAndEnergy( numParticles, allPosData, allParticleIndices, *neighborList, forceData);

    return static_cast<double>(energy);
}

void ReferenceCalcMBPolThreeBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolThreeBodyForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {

        std::vector<int> particleIndices;
        force.getParticleParameters(i, particleIndices);
        allParticleIndices[i] = particleIndices;

    }
}

ReferenceCalcMBPolDispersionForceKernel::ReferenceCalcMBPolDispersionForceKernel(std::string name, const Platform& platform, const OpenMM::System& system) :
       CalcMBPolDispersionForceKernel(name, platform), system(system) {
    useCutoff = 0;
    usePBC = 0;
    cutoff = 1.0e+10;
    neighborList = NULL;
}

ReferenceCalcMBPolDispersionForceKernel::~ReferenceCalcMBPolDispersionForceKernel() {
    if( neighborList ){
        delete neighborList;
    }
}

void ReferenceCalcMBPolDispersionForceKernel::initialize(const OpenMM::System& system, const MBPolDispersionForce& force) {

    // per-particle parameters

    numParticles = force.getNumMolecules();
    allParticleIndices.resize(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){

        std::vector<int> particleIndices;
        force.getParticleParameters(ii, particleIndices );
        allParticleIndices[ii] = particleIndices;

    }

    useCutoff              = (force.getNonbondedMethod() != MBPolDispersionForce::NoCutoff);
    usePBC                 = (force.getNonbondedMethod() == MBPolDispersionForce::CutoffPeriodic);
    cutoff                 = force.getCutoff();
    neighborList           = useCutoff ? new NeighborList() : NULL;

}

double ReferenceCalcMBPolDispersionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& allPosData   = extractPositions(context);
    vector<RealVec> posData;
    posData.resize(numParticles);
    vector<set<int> > allExclusions;
    allExclusions.resize(numParticles);
    // posData has only oxygens
    for( int ii = 0; ii < numParticles; ii++ ){
        posData[ii] = allPosData[allParticleIndices[ii][0]];
    }
    vector<RealVec>& forceData = extractForces(context);
    MBPolReferenceDispersionForce dispersionForce;
    RealOpenMM energy;
    dispersionForce.setCutoff( cutoff );
    // neighborList created only with oxygens, then allParticleIndices is used to get reference to the hydrogens
    computeNeighborListVoxelHash( *neighborList, numParticles, posData, allExclusions, extractBoxSize(context), usePBC, cutoff, 0.0, false);
    if( usePBC ){
        dispersionForce.setNonbondedMethod( MBPolReferenceDispersionForce::CutoffPeriodic);
        RealVec& box = extractBoxSize(context);
        double minAllowedSize = 1.999999*cutoff;
        if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
            throw OpenMMException("The periodic box size has decreased to less than twice the cutoff.");
        }
        dispersionForce.setPeriodicBox(box);
    } else {
        dispersionForce.setNonbondedMethod( MBPolReferenceDispersionForce::CutoffNonPeriodic);
    }
    // here we need allPosData, every atom!
    energy  = dispersionForce.calculateForceAndEnergy( numParticles, allPosData, allParticleIndices, *neighborList, forceData);

    return static_cast<double>(energy);
}

void ReferenceCalcMBPolDispersionForceKernel::copyParametersToContext(ContextImpl& context, const MBPolDispersionForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {

        std::vector<int> particleIndices;
        force.getParticleParameters(i, particleIndices);
        allParticleIndices[i] = particleIndices;

    }
}
