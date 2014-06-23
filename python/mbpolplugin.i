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

namespace MBPolPlugin {

class MBPolElectrostaticsForce : public OpenMM::Force {
public:
    MBPolElectrostaticsForce();
    int getNumElectrostaticss() const;
    double getCutoffDistance() const;
    void setCutoffDistance(double distance);
    void updateParametersInContext(OpenMM::Context& context);
};

}

