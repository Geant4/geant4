// $Id: G4ChordFinderSaf.hh,v 1.1 2003-11-13 19:05:26 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4ChordFinderRad
//
// Class description:
//
// A class that provides the next chord, utilising an estimate of
//  the radius of curvature.

// History:
// - 07.11.03 John Apostolakis,  design and implementation 
// -------------------------------------------------------------------

#include "G4ChordFinder.hh"

class G4ChordFinderSaf : public G4ChordFinder
{

public:
  G4ChordFinderSaf(G4MagInt_Driver* pIntegrationDriver); 
    
  G4ChordFinderSaf( G4MagneticField*        theMagField,
		    G4double                stepMinimum, 
		    G4MagIntegratorStepper* pItsStepper ); 

  ~G4ChordFinderSaf(); 

  G4double FindNextChord( const  G4FieldTrack  yStart,
			  G4double     stepMax,
			  G4FieldTrack&   yEnd,  // Endpoint
			  G4double&   dyErrPos,  // Error of endpoint
			  G4double    epsStep,
			  G4double*  pStepForAccuracy,
			  const G4ThreeVector latestSafetyOrigin,
			  G4double       latestSafetyRadius 
			  );  

  void PrintStatistics();

private:
  // G4int fNoInitialRadBig,  fNoInitialRadSmall; 
  // G4int fNoTrialsRadBig,   fNoTrialsRadSmall; 

};
