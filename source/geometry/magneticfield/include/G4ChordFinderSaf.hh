//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4ChordFinderSaf.hh,v 1.2 2003/12/09 15:35:07 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
