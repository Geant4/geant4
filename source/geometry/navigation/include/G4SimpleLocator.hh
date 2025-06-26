//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4SimpleLocator 
//
// Class description:
// 
// Implementing the calculation of the intersection point with a boundary when
// PropagationInField is used. Derived from method LocateIntersectionPoint()
// from G4PropagatorInField, it is based on a linear method for finding the
// intersection point; the difference compared to G4MultiLevelLocator is that
// no 'depth' algorithm is used in case of slow progress for finding the
// intersection point.

// Author: Tatiana Nikitina (CERN), 27 October 20008.
// ---------------------------------------------------------------------------
#ifndef G4SIMPLELOCATOR_HH
#define G4SIMPLELOCATOR_HH 1

#include "G4VIntersectionLocator.hh"

/**
 * @brief G4SimpleLocator implements the calculation of the intersection point
 * with a boundary when G4PropagationInField is used. It is based on a linear
 * method for finding the intersection point; the difference compared to
 * G4MultiLevelLocator is that no 'depth' algorithm is used in case of slow
 * progress for finding the intersection point.
 */

class G4SimpleLocator : public G4VIntersectionLocator
{
   public:
 
    /**
     * Constructor and default Destructor.
     */
     G4SimpleLocator(G4Navigator* aNavigator);
     ~G4SimpleLocator() override;

    /**
     * If such an intersection exists, this function calculates the
     * intersection point of the true path of the particle with the surface
     * of the current volume (or of one of its daughters). 
     *  @note Should use lateral displacement as measure of convergence.
     *  @param[in] curveStartPointTangent Start point tangent track.
     *  @param[in] curveEndPointTangent End point tangent track.
     *  @param[in] trialPoint Trial point.
     *  @param[out] intersectPointTangent Intersection point tangent track.
     *  @param[out] recalculatedEndPoint Flagging if end point was recomputed.
     *  @param[in,out] fPreviousSafety Previous safety distance.
     *  @param[in,out] fPreviousSftOrigin Previous safety point origin.
     *  @returns Whether intersection exists or not. 
     */
     G4bool EstimateIntersectionPoint( 
         const  G4FieldTrack&       curveStartPointTangent,           // A
         const  G4FieldTrack&       curveEndPointTangent,             // B
         const  G4ThreeVector&      trialPoint,                       // E
                G4FieldTrack&       intersectPointTangent,            // Output
	        G4bool&             recalculatedEndPoint,             // Out
                G4double&           fPreviousSafety,                  // In/Out
                G4ThreeVector&      fPreviousSftOrigin) override;     // In/Out
};

#endif
