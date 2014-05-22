/* -------------------------------------------------------------------------- *
 *                                OpenMMMBPol                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/MBPolOneBodyForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerMBPolSerializationProxies();

void testSerialization() {
    // Create a Force.

    MBPolOneBodyForce force1;
    force1.addOneBody(0, 1, 3);
    force1.addOneBody(2, 4, 4);
    force1.addOneBody(5, 0, 1);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<MBPolOneBodyForce>(&force1, "Force", buffer);
    MBPolOneBodyForce* copy = XmlSerializer::deserialize<MBPolOneBodyForce>(buffer);

    // Compare the two forces to see if they are identical.  
    MBPolOneBodyForce& force2 = *copy;
    ASSERT_EQUAL(force1.getNumOneBodys(), force2.getNumOneBodys());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumOneBodys()); ii++) {
        int p11, p12, p13;
        int p21, p22, p23;

        force1.getOneBodyParameters(ii, p11, p12, p13);
        force2.getOneBodyParameters(ii, p21, p22, p23);

        ASSERT_EQUAL(p11, p21);
        ASSERT_EQUAL(p12, p22);
        ASSERT_EQUAL(p13, p23);
    }
}

int main() {
    try {
        registerMBPolSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

