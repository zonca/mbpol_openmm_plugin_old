
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "MBPolReferenceForce.h"
#include <vector>

using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   Load delta of two vectors

   @param xVector      first vector
   @param yVector      second vector
   @param deltaR       output vector: y - x

   --------------------------------------------------------------------------------------- */

void MBPolReferenceForce::loadDeltaR( const RealVec& xVector, const RealVec& yVector,
                                       std::vector<RealOpenMM>& deltaR ){


   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::loadDeltaR";

   // ---------------------------------------------------------------------------------------

    deltaR.resize(0);
    deltaR.push_back( yVector[0] - xVector[0] );
    deltaR.push_back( yVector[1] - xVector[1] );
    deltaR.push_back( yVector[2] - xVector[2] );
}

/**---------------------------------------------------------------------------------------

   Calculate norm squared of 3d vector

   @param inputVector    vector whose norm squared is to be computed

   @return norm squared

   --------------------------------------------------------------------------------------- */

RealOpenMM MBPolReferenceForce::getNormSquared3( const std::vector<RealOpenMM>& inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getNorm3";

   // ---------------------------------------------------------------------------------------

   // get 3 norm

   return ( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
}

/**---------------------------------------------------------------------------------------

   Calculate norm squared of 3d vector

   @param inputVector    vector whose norm squared is to be computed

   @return norm squared

   --------------------------------------------------------------------------------------- */

RealOpenMM MBPolReferenceForce::getNormSquared3( const RealOpenMM* inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getNorm3";

   // ---------------------------------------------------------------------------------------

   // get 3 norm

   return ( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
}

/**---------------------------------------------------------------------------------------

   Calculate norm of 3d vector

   @param inputVector            vector whose norm is to be computed

   @return norm

   --------------------------------------------------------------------------------------- */

RealOpenMM MBPolReferenceForce::getNorm3( const std::vector<RealOpenMM>& inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getNorm3";

   // ---------------------------------------------------------------------------------------

   // get 3 norm

   return SQRT( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
}

RealOpenMM MBPolReferenceForce::getNorm3( const RealOpenMM* inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getNorm3";

   // ---------------------------------------------------------------------------------------

   // get 3 norm

   return SQRT( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
}

RealOpenMM MBPolReferenceForce::normalizeVector3( RealOpenMM* inputVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::normalizeVector3";

   // ---------------------------------------------------------------------------------------

    RealOpenMM norm   = SQRT( inputVector[0]*inputVector[0] + inputVector[1]*inputVector[1] + inputVector[2]*inputVector[2] );
    if( norm > 0.0 ){
        RealOpenMM normI  = 1.0/norm;
        inputVector[0]   *= normI;
        inputVector[1]   *= normI;
        inputVector[2]   *= normI;
    }

    return norm;
}

/**---------------------------------------------------------------------------------------

   Calculate dot product of 3d vectors

   @param xVector   first vector
   @param yVector   second vector

   @return dot product

   --------------------------------------------------------------------------------------- */

RealOpenMM MBPolReferenceForce::getDotProduct3( const std::vector<RealOpenMM>& xVector, const std::vector<RealOpenMM>& yVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getDotProduct3";

   // ---------------------------------------------------------------------------------------

   // get dot product

   return xVector[0]*yVector[0] + xVector[1]*yVector[1] + xVector[2]*yVector[2];
}

/**---------------------------------------------------------------------------------------

   Calculate dot product of 3d vectors

   @param xVector   first vector
   @param yVector   second vector

   @return dot product

   --------------------------------------------------------------------------------------- */

RealOpenMM MBPolReferenceForce::getDotProduct3( const RealOpenMM* xVector, const RealOpenMM* yVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getDotProduct3";

   // ---------------------------------------------------------------------------------------

   // get dot product

   return xVector[0]*yVector[0] + xVector[1]*yVector[1] + xVector[2]*yVector[2];
}

RealOpenMM MBPolReferenceForce::getDotProduct3( const RealOpenMM* xVector, const OpenMM::Vec3& yVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getDotProduct3";

   // ---------------------------------------------------------------------------------------

   // get dot product

   return xVector[0]*yVector[0] + xVector[1]*yVector[1] + xVector[2]*yVector[2];
}

/**---------------------------------------------------------------------------------------

   Calculate dot product of 3d vectors

   @param vectorOffset offset into first first vector
   @param xVector      first vector
   @param yVector      second vector

   @return dot product

   --------------------------------------------------------------------------------------- */

RealOpenMM MBPolReferenceForce::getDotProduct3( unsigned int vectorOffset, const std::vector<RealOpenMM>& xVector, const RealOpenMM* yVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getDotProduct3";

   // ---------------------------------------------------------------------------------------

   // get dot product

   return xVector[vectorOffset+0]*yVector[0] + xVector[vectorOffset+1]*yVector[1] + xVector[vectorOffset+2]*yVector[2];
}

/**---------------------------------------------------------------------------------------

   Calculate z = x X y

   @param xVector      input vector
   @param yVector      input vector
   @param zVector      output vector: z = x X y

   --------------------------------------------------------------------------------------- */

void MBPolReferenceForce::getCrossProduct( const std::vector<RealOpenMM>& xVector,
                                            const std::vector<RealOpenMM>& yVector,
                                            std::vector<RealOpenMM>& zVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getCrossProduct";

   // ---------------------------------------------------------------------------------------

   zVector[0]  = xVector[1]*yVector[2] - xVector[2]*yVector[1];
   zVector[1]  = xVector[2]*yVector[0] - xVector[0]*yVector[2];
   zVector[2]  = xVector[0]*yVector[1] - xVector[1]*yVector[0];

   return;
}

/**---------------------------------------------------------------------------------------

   Calculate z = x X y

   @param xVector      input vector
   @param yVector      input vector
   @param zVector      output vector: z = x X y

   --------------------------------------------------------------------------------------- */

void MBPolReferenceForce::getCrossProduct( const RealOpenMM* xVector,
                                            const RealOpenMM* yVector,
                                            RealOpenMM* zVector ){

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "MBPolReferenceForce::getCrossProduct";

   // ---------------------------------------------------------------------------------------

   zVector[0]  = xVector[1]*yVector[2] - xVector[2]*yVector[1];
   zVector[1]  = xVector[2]*yVector[0] - xVector[0]*yVector[2];
   zVector[2]  = xVector[0]*yVector[1] - xVector[1]*yVector[0];

   return;
}

