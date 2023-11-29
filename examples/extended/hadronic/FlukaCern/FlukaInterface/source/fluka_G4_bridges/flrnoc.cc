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
#ifdef G4_USE_FLUKA


#include "flrnoc.hh"

#include "Randomize.hh"
#include "G4Exception.hh"


void flrnoc_(const G4int& /*dummySeed1*/, const G4int& /*dummySeed2*/, 
	     G4int& randomCallMod1BCounter, const G4int& /*dummyRandomCallBillionCounter*/) {

	randomCallMod1BCounter = 1;

	G4Exception("flrnoc_",
		    "Calls to flrnoc are not supported.\n" \
		    "Would require access to G4 random engine internal status (number of generated numbers).\n" \
		    "This would be overkill, because:\n" \
		    "(1) Display of random engine internal status should be done through FL64WR anyway (calling G4Random).\n" \
		    "(2) External use of random engine internal status is bad practice. To trace back where a call is from, could just use booleans passed as arguments.\n",
		    FatalException,
		    "Unsupported function call.");

}


#endif // G4_USE_FLUKA
