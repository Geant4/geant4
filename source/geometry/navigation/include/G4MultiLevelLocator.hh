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
// Class description:
// 
// Calculation of Intersection Point with a boundary 
// when PropagationInField is used.
// Derived from method LocateIntersectionPoint(...) from G4PropagatorInField
// It is based on Linear Method for Finding Intersection Point with
// a 'depth' algorithm in case of slow progress.(Intersection is not found
// after 10 trials)
//
//
// History:
// -------
// 27.10.08  Tatiana Nikitina  derived from LocateIntersectionPoint(...) from 
//                             G4PropagatorInField class
//
#ifndef G4MULTILEVELLOCATOR_HH
#define G4MULTILEVELLOCATOR_HH

#include "G4VIntersectionLocator.hh"

class G4MultiLevelLocator:public G4VIntersectionLocator
 {
   public:  // with description 
 
     G4MultiLevelLocator(G4Navigator *theNavigator);
        // Constructor
     ~G4MultiLevelLocator();
        // Default destructor.
     
     G4bool EstimateIntersectionPoint( 
         const  G4FieldTrack&       curveStartPointTangent,  //  A
         const  G4FieldTrack&       curveEndPointTangent,    //  B
         const  G4ThreeVector&      trialPoint,              //  E
                G4FieldTrack&       intersectPointTangent,   // Output
	        G4bool&             recalculatedEndPoint,    // Out
                G4double&           fPreviousSafety=0,//In/Out
                G4ThreeVector&      fPreviousSftOrigin=0);//In/Out

                                                        
      // If such an intersection exists, this function 
      // calculate the intersection point of the true path of the particle 
      // with the surface of the current volume (or of one of its daughters). 
      // (Should use lateral displacement as measure of convergence). 

    private:

   static const G4int max_depth=4;
   // static const G4int max_depth=7;
   G4FieldTrack* ptrInterMedFT[max_depth+1];
     // Used to store intermediate values of tracks in case of
     // too slow progress

   G4int maxNumberOfStepsForIntersection;
   G4int maxNumberOfCallsToReIntegration;
   G4int maxNumberOfCallsToReIntegration_depth;
    //  Counters for Statistics about Location and ReIntegrations
  

 };

#endif
