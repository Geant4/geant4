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
//
// $Id: G4SimpleLocator.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// Class G4SimpleLocator 
//
// class description:
// 
// Implementing the calculation of the intersection point with a boundary when
// PropagationInField is used. Derived from method LocateIntersectionPoint()
// from G4PropagatorInField, it is based on a linear method for finding the
// intersection point; the difference compared to G4MultiLevelLocator is that
// no 'depth' algorithm is used in case of slow progress for finding the
// intersection point.

// History:
// -------
// 27.10.08 - Tatiana Nikitina: Derived from LocateIntersectionPoint() from 
//                              G4PropagatorInField class
// ---------------------------------------------------------------------------

#ifndef G4SIMPLELOCATOR_HH
#define G4SIMPLELOCATOR_HH

#include "G4VIntersectionLocator.hh"

class G4SimpleLocator : public G4VIntersectionLocator
{
   public:  // with description 
 
     G4SimpleLocator(G4Navigator *theNavigator);
        // Constructor
     ~G4SimpleLocator();
        // Default destructor

     G4bool EstimateIntersectionPoint( 
         const  G4FieldTrack&       curveStartPointTangent,  // A
         const  G4FieldTrack&       curveEndPointTangent,    // B
         const  G4ThreeVector&      trialPoint,              // E
                G4FieldTrack&       intersectPointTangent,   // Output
	        G4bool&             recalculatedEndPoint,    // Out
                G4double&           fPreviousSafety,         // In/Out
                G4ThreeVector&      fPreviousSftOrigin);     // In/Out
      // If such an intersection exists, this function calculates the
      // intersection point of the true path of the particle with the surface
      // of the current volume (or of one of its daughters). 
      // Should use lateral displacement as measure of convergence
};

#endif
