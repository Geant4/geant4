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
// $Id: G4ChordFinderSaf.hh 107059 2017-11-01 14:58:16Z gcosmo $
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
  G4ChordFinderSaf(G4VIntegrationDriver* pIntegrationDriver);
    
  G4ChordFinderSaf( G4MagneticField*        theMagField,
		    G4double                stepMinimum, 
		    G4MagIntegratorStepper* pItsStepper ); 

  ~G4ChordFinderSaf(); 

  G4double FindNextChord( const  G4FieldTrack&  yStart,
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
