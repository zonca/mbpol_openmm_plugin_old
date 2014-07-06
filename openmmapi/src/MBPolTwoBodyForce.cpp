/* -------------------------------------------------------------------------- *
 *                                OpenMMMBPol                                *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/MBPolTwoBodyForce.h"
#include "openmm/internal/MBPolTwoBodyForceImpl.h"
#include <iostream>

using namespace  OpenMM;
using namespace MBPolPlugin;
using std::string;
using std::vector;

MBPolTwoBodyForce::MBPolTwoBodyForce() : nonbondedMethod(CutoffNonPeriodic), cutoff(1.0e+10) {
}

int MBPolTwoBodyForce::addParticle(const std::vector<int> & particleIndices ) {
    std::cout << particleIndices[0] << std::endl;
    parameters.push_back(TwoBodyInfo(particleIndices));
    return parameters.size()-1;
}

int MBPolTwoBodyForce::getNumMolecules() const {
    return parameters.size();
}

void MBPolTwoBodyForce::getParticleParameters(int particleIndex, std::vector<int>& particleIndices ) const {
    particleIndices     = parameters[particleIndex].particleIndices;
}

void MBPolTwoBodyForce::setParticleParameters(int particleIndex, std::vector<int>& particleIndices  ) {
      parameters[particleIndex].particleIndices =particleIndices;
}

void MBPolTwoBodyForce::setCutoff( double inputCutoff ){
    cutoff = inputCutoff;
}

double MBPolTwoBodyForce::getCutoff( void ) const {
    return cutoff;
}

MBPolTwoBodyForce::NonbondedMethod MBPolTwoBodyForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void MBPolTwoBodyForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

ForceImpl* MBPolTwoBodyForce::createImpl() const {
    return new MBPolTwoBodyForceImpl(*this);
}

void MBPolTwoBodyForce::updateParametersInContext(Context& context) {
    dynamic_cast<MBPolTwoBodyForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
