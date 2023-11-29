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
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDummySD  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************
//
// Dummy sensitive used only to flag sensitivity
// in cells of RO geometry.
//

#ifndef GammaRayTelDummySD_h
#define GammaRayTelDummySD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;

class [[deprecated]] GammaRayTelDummySD: public G4VSensitiveDetector {
public:
    explicit GammaRayTelDummySD(G4String name) : G4VSensitiveDetector(name) {
	}

	~GammaRayTelDummySD() override;

	void Initialize(G4HCofThisEvent *collection);

	auto ProcessHits(G4Step *step, G4TouchableHistory *history) -> G4bool {
		return false;
	}

	void EndOfEvent(G4HCofThisEvent *collection) {
	}

	void clear() {
	}

	void DrawAll() {
	}

	void PrintAll() {
	}
};
#endif
