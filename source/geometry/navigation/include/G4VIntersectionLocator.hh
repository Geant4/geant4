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
// Gives possibility to choose the method of Intersection :
// G4SimpleLocator, G4MultiLevelLocator, G4BrentLocator
//
// Key Method:
//              EstimateIntersectionPoint(..)
// History:
// -------
// 27.10.07 John Apostolakis, Tatiana Nikitina  design and implementation 
//

#ifndef G4VINTERSECTIONLOCATOR_HH
#define G4VINTERSECTIONLOCATOR_HH

#include "G4Types.hh" 
#include "G4ThreeVector.hh"
#include "G4FieldTrack.hh"

#include "G4Navigator.hh"
#include "G4ChordFinder.hh"

class G4VIntersectionLocator
 {
   public:  // with description 
 
     G4VIntersectionLocator(G4Navigator *theNavigator);
        // Constructor
     virtual ~G4VIntersectionLocator();
        // Default destructor.
     
     virtual G4bool EstimateIntersectionPoint( 
         const  G4FieldTrack&       curveStartPointTangent,  //  A
         const  G4FieldTrack&       curveEndPointTangent,    //  B
         const  G4ThreeVector&      trialPoint,              //  E
                G4FieldTrack&       intersectPointTangent,   // Output
	        G4bool&             recalculatedEndPoint,    // Out
                G4double&           fPreviousSafety,//In/Out
                G4ThreeVector&      fPreviousSftOrigin)//In/Out)
                                                        = 0;    
      // If such an intersection exists, this function 
      // calculate the intersection point of the true path of the particle 
      // with the surface of the current volume (or of one of its daughters). 
      // (Should use lateral displacement as measure of convergence). 
      // (NOTE: Changes Safety)
     
     void printStatus( const G4FieldTrack&        startFT,
                       const G4FieldTrack&        currentFT, 
                             G4double             requestStep, 
                             G4double             safety,
                             G4int                step); 
    // Print Method - useful mostly for debugging.
     
  public:  // without dascription
   // inline Set/Get methods used for IntersctionLocator
   inline G4double  GetDeltaIntersectionFor(){return fiDeltaIntersection;} 
   inline G4double  GetEpsilonStepFor() {return fiEpsilonStep;}
   inline G4Navigator*  GetNavigatorFor(){return fiNavigator;} 
   inline G4ChordFinder*  GetChordFinderFor() {return fiChordFinder;}
   inline G4int  GetVerboseFor() {return fVerboseLevel;}

   inline void SetEpsilonStepFor( G4double EpsilonStep ){fiEpsilonStep=EpsilonStep;}  
   inline void SetDeltaIntersectionFor( G4double deltaIntersection ){
                                        fiDeltaIntersection=deltaIntersection;}
   inline void SetNavigatorFor( G4Navigator *fNavigator ){fiNavigator=fNavigator;}  
   inline void SetChordFinderFor(G4ChordFinder *fCFinder ){fiChordFinder=fCFinder;}
   inline void SetSafetyParametersFor(G4bool UseSafety ){fiUseSafety=UseSafety;}
   inline void SetVerboseFor(G4int fVerbose){fVerboseLevel=fVerbose;}

   //with description
   inline G4bool IntersectChord( G4ThreeVector  StartPointA, 
                                 G4ThreeVector  EndPointB,
                                 G4double      &NewSafety,
                                 G4double      &fPreviousSafety,//In/Out
		                 G4ThreeVector &fPreviousSftOrigin,//In/Out
                                 G4double      &LinearStepLength,
                                 G4ThreeVector &IntersectionPoint);
      // Intersect the chord from StartPointA to EndPointB
      // and return whether an intersection occurred
      // (NOTE : changes Safety)

 protected:  // with description

 G4FieldTrack ReEstimateEndpoint( const G4FieldTrack &CurrentStateA,  
                                    const G4FieldTrack &EstimtdEndStateB,
                                          G4double      linearDistSq,
                                          G4double      curveDist);
     // Return new estimate for state after curveDist 
     // starting from CurrentStateA,  to replace EstimtdEndStateB,
     // (and report displacement -- if field is compiled verbose.)
  
 protected:

   G4double kCarTolerance;
   G4int    fVerboseLevel;
   // For verbose purposes

   G4double      fiEpsilonStep;
   G4double      fiDeltaIntersection;
   G4Navigator   *fiNavigator;
   G4ChordFinder *fiChordFinder;
   G4bool        fiUseSafety;
   // For passing the parameters from G4PropagatorInField
  
 };

inline G4bool G4VIntersectionLocator::IntersectChord(
                          G4ThreeVector  StartPointA, 
                          G4ThreeVector  EndPointB,
                          G4double      &NewSafety,
                          G4double      &fPreviousSafety,//In/Out
		          G4ThreeVector &fPreviousSftOrigin,//In/Out
                          G4double      &LinearStepLength,
                          G4ThreeVector &IntersectionPoint)

{
    // Calculate the direction and length of the chord AB
    G4ThreeVector  ChordAB_Vector = EndPointB - StartPointA;
    G4double       ChordAB_Length = ChordAB_Vector.mag();  // Magnitude (norm)
    G4ThreeVector  ChordAB_Dir =    ChordAB_Vector.unit();
    G4bool intersects;
    G4ThreeVector OriginShift = StartPointA -  fPreviousSftOrigin ;
    G4double      MagSqShift  = OriginShift.mag2() ;
    G4double      currentSafety;
    G4bool        doCallNav= false;

    if( MagSqShift >= sqr(fPreviousSafety) )
    {
        currentSafety = 0.0 ;
    }else{
        currentSafety = fPreviousSafety - std::sqrt(MagSqShift) ;
    }

    if( fiUseSafety && (ChordAB_Length <= currentSafety) )
    {
       // The Step is guaranteed to be taken
       LinearStepLength = ChordAB_Length;
       intersects = false;
       NewSafety= currentSafety;
    }
    else
    {
       doCallNav= true; 
       // Check whether any volumes are encountered by the chord AB

       // G4cout << " G4PropagatorInField calling Navigator::ComputeStep " << G4endl ;
       LinearStepLength = 
        GetNavigatorFor()->ComputeStep( StartPointA, ChordAB_Dir,
                                 ChordAB_Length, NewSafety );
       intersects = (LinearStepLength <= ChordAB_Length); 
       // G4Navigator contracts to return k_infinity if len==asked
       // and it did not find a surface boundary at that length
       LinearStepLength = std::min( LinearStepLength, ChordAB_Length);

       // G4cout << " G4PiF got step= " << LinearStepLength << " safety= " << NewSafety << G4endl;

       // Save the last calculated safety!
       fPreviousSftOrigin = StartPointA;
       fPreviousSafety= NewSafety;

       if( intersects ){
          // Intersection Point of chord AB and either volume A's surface 
          //                                or a daughter volume's surface ..
          IntersectionPoint = StartPointA + LinearStepLength * ChordAB_Dir;
       }
    }
 
    return intersects;
}

#endif
