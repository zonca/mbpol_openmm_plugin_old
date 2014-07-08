/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/MBPolDispersionForceImpl.h"
#include "openmm/mbpolKernels.h"
#include <map>
#include <cmath>

using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;

using std::pair;
using std::vector;
using std::set;

MBPolDispersionForceImpl::MBPolDispersionForceImpl(const MBPolDispersionForce& owner) : owner(owner) {
}

MBPolDispersionForceImpl::~MBPolDispersionForceImpl() {
}

void MBPolDispersionForceImpl::initialize(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();

    // check that cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == MBPolDispersionForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoff();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("MBPolDispersionForce: The cutoff distance cannot be greater than half the periodic box size.");
    }   

    kernel = context.getPlatform().createKernel(CalcMBPolDispersionForceKernel::Name(), context);
    kernel.getAs<CalcMBPolDispersionForceKernel>().initialize(context.getSystem(), owner);
}

double MBPolDispersionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcMBPolDispersionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> MBPolDispersionForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcMBPolDispersionForceKernel::Name());
    return names;
}

void MBPolDispersionForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcMBPolDispersionForceKernel>().copyParametersToContext(context, owner);
}
