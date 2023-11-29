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
//      ------------ GammaRayTelPrimaryGeneratorAction  ------
//           by G.Santin, F.Longo & R.Giannitrapani (30 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelPrimaryGeneratorAction_h
#define GammaRayTelPrimaryGeneratorAction_h 1

#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class GammaRayTelDetectorConstruction;
class GammaRayTelPrimaryGeneratorMessenger;
class G4GeneralParticleSource;

class GammaRayTelPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
	explicit GammaRayTelPrimaryGeneratorAction();

	~GammaRayTelPrimaryGeneratorAction() override;

	void GeneratePrimaries(G4Event *event) override;

	void SetRndmFlag(G4String value) {
		rndmFlag = value;
	}

	void SetSourceType(G4int value) {
		sourceType = value;
	}

	void SetSpectrumType(G4int value) {
		spectrumType = value;
	}

	void SetVertexRadius(G4double value) {
		vertexRadius = value;
	}

	void SetSourceGen(G4bool value) {
		sourceGun = value;
	}

private:
	G4ParticleGun *particleGun;

	G4GeneralParticleSource *particleSource;

	const GammaRayTelDetectorConstruction *detector;

	GammaRayTelPrimaryGeneratorMessenger *gunMessenger;

	G4String rndmFlag{"off"}; // flag for a random impact point

	G4int sourceType{0};

	G4double vertexRadius{15. * cm};

	G4int spectrumType{0};

	G4bool sourceGun{false}; // false for GeneralParticleSource
};
#endif
