%module mbpolplugin

%import(module="simtk.openmm") "OpenMMSwigHeaders.i"

%{
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/MBPolOneBodyForce.h"
#include "openmm/MBPolTwoBodyForce.h"
#include "openmm/MBPolThreeBodyForce.h"
#include "openmm/MBPolDispersionForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
%}

// %feature("autodoc", "1");
// %nodefaultctor;

using namespace OpenMM;

namespace MBPolPlugin {

class MBPolElectrostaticsForce : public OpenMM::Force {
public:
    MBPolElectrostaticsForce();

    int getNumElectrostatics() const;

    // NonbondedMethod getNonbondedMethod() const;

    // void setNonbondedMethod(NonbondedMethod method);

    // PolarizationType getPolarizationType() const;

    // void setPolarizationType(PolarizationType type);

    // double getCutoffDistance(void) const;

    // void setCutoffDistance(double distance);

    // void setIncludeChargeRedistribution( bool chargeRedistribution );

    // bool getIncludeChargeRedistribution( void ) const;

    // double getAEwald() const;

    // void setAEwald(double aewald);

    // int getPmeBSplineOrder() const;

    // void getPmeGridDimensions(std::vector<int>& gridDimension) const;

    // void setPmeGridDimensions(const std::vector<int>& gridDimension);

    // int addElectrostatics(double charge, const std::vector<double>& molecularDipole, const std::vector<double>& molecularQuadrupole, int axisType,
    //                  int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, const std::vector<double>& thole, double dampingFactor, double polarity);

    // void getElectrostaticsParameters(int index, double& charge, std::vector<double>& molecularDipole, std::vector<double>& molecularQuadrupole,
    //                             int& axisType, int& multipoleAtomZ, int& multipoleAtomX, int& multipoleAtomY, std::vector<double>& thole, double& dampingFactor, double& polarity) const;

    // void setElectrostaticsParameters(int index, double charge, const std::vector<double>& molecularDipole, const std::vector<double>& molecularQuadrupole,
    //                             int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, const std::vector<double>& thole, double dampingFactor, double polarity);

    // void setCovalentMap(int index, CovalentType typeId, const std::vector<int>& covalentAtoms);

    // void getCovalentMap(int index, CovalentType typeId, std::vector<int>& covalentAtoms) const;

    // void getCovalentMaps(int index, std::vector < std::vector<int> >& covalentLists) const;

    // int getMutualInducedMaxIterations(void) const;

    // void setMutualInducedMaxIterations(int inputMutualInducedMaxIterations);

    // double getMutualInducedTargetEpsilon(void) const;

    // void setMutualInducedTargetEpsilon(double inputMutualInducedTargetEpsilon);

    // double getEwaldErrorTolerance() const;

    // void setEwaldErrorTolerance(double tol);

    // void getElectrostaticPotential(const std::vector< Vec3 >& inputGrid,
    //                                 Context& context, std::vector< double >& outputElectrostaticPotential);

    // void updateParametersInContext(Context& context);
};

}
