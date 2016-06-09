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
// $Id: G4VIntersectionLocator.hh,v 1.4 2008/11/14 18:26:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// Class G4VIntersectionLocator 
//
// class description:
// 
// Base class for the calculation of the intersection point with a boundary 
// when PropagationInField is used.
// Gives possibility to choose the method of intersection; concrete locators
// implemented are: G4SimpleLocator, G4MultiLevelLocator, G4BrentLocator.
//
// Key Method: EstimateIntersectionPoint()

// History:
// -------
// 27.10.08 - John Apostolakis, Tatiana Nikitina: Design and implementation 
// ---------------------------------------------------------------------------

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
       // Default destructor
     
     virtual G4bool EstimateIntersectionPoint( 
         const  G4FieldTrack&       curveStartPointTangent,  // A
         const  G4FieldTrack&       curveEndPointTangent,    // B
         const  G4ThreeVector&      trialPoint,              // E
                G4FieldTrack&       intersectPointTangent,   // Output
                G4bool&             recalculatedEndPoint,    // Out
                G4double&           fPreviousSafety,         // In/Out
                G4ThreeVector&      fPreviousSftOrigin) = 0; // In/Out   
       // If such an intersection exists, this function calculates the
       // intersection point of the true path of the particle with the surface
       // of the current volume (or of one of its daughters). 
       // Should use lateral displacement as measure of convergence
       // NOTE: changes the safety!

     void printStatus( const G4FieldTrack& startFT,
                       const G4FieldTrack& currentFT, 
                             G4double      requestStep, 
                             G4double      safety,
                             G4int         step);
       // Print Method, useful mostly for debugging

     inline G4bool IntersectChord( G4ThreeVector  StartPointA, 
                                   G4ThreeVector  EndPointB,
                                   G4double      &NewSafety,
                                   G4double      &fPreviousSafety,    // In/Out
                                   G4ThreeVector &fPreviousSftOrigin, // In/Out
                                   G4double      &LinearStepLength,
                                   G4ThreeVector &IntersectionPoint);
       // Intersect the chord from StartPointA to EndPointB and return
       // whether an intersection occurred. NOTE: changes the Safety!

  public:  // without description

    // inline Set/Get methods used for IntersctionLocator

    inline G4double  GetDeltaIntersectionFor();
    inline G4double  GetEpsilonStepFor();
    inline G4Navigator*  GetNavigatorFor();
    inline G4ChordFinder*  GetChordFinderFor();
    inline G4int  GetVerboseFor();

    inline void SetEpsilonStepFor( G4double EpsilonStep );
    inline void SetDeltaIntersectionFor( G4double deltaIntersection );
    inline void SetNavigatorFor( G4Navigator *fNavigator );
    inline void SetChordFinderFor(G4ChordFinder *fCFinder );
    inline void SetSafetyParametersFor(G4bool UseSafety );
    inline void SetVerboseFor(G4int fVerbose);

    inline void AddAdjustementOfFoundIntersection(G4bool UseCorrection);
    inline G4bool GetAdjustementOfFoundIntersection();

  protected:  // with description

    G4FieldTrack ReEstimateEndpoint( const G4FieldTrack &CurrentStateA,  
                                     const G4FieldTrack &EstimtdEndStateB,
                                           G4double      linearDistSq,
                                           G4double      curveDist);
      // Return new estimate for state after curveDist starting from
      // CurrentStateA, to replace EstimtdEndStateB, and report displacement
      // (if field is compiled verbose)

     
    G4ThreeVector GetLocalSurfaceNormal(const G4ThreeVector &CurrentE_Point,
                                              G4bool &validNormal);
      // Return the SurfaceNormal of Intersecting Solid

    G4bool AdjustmentOfFoundIntersection(const G4ThreeVector &A,
                                         const G4ThreeVector &CurrentE_Point, 
                                         const G4ThreeVector &CurrentF_Point,
                                         const G4ThreeVector &MomentumDir,
                                         const G4bool         IntersectAF, 
                                               G4ThreeVector &IntersectionPoint,
                                               G4double      &NewSafety,
                                               G4double      &fPrevSafety,
                                               G4ThreeVector &fPrevSftOrigin );
      // Optional method for adjustment of located intersection point
      // using the surface-normal
  
  protected:

    G4double kCarTolerance;
    G4int    fVerboseLevel;
      // For verbose purposes
    G4bool   fUseNormalCorrection;

    G4double       fiEpsilonStep;
    G4double       fiDeltaIntersection;
    G4Navigator   *fiNavigator;
    G4ChordFinder *fiChordFinder;
    G4bool         fiUseSafety;
      // For passing the parameters from G4PropagatorInField

    G4Navigator *fHelpingNavigator;
      // Helper for location
};

#include "G4VIntersectionLocator.icc"

#endif
