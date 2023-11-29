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
//      ------------ GammaRayTelCalorimeterSD  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelCalorimeterSD_h
#define GammaRayTelCalorimeterSD_h 1

#include "GammaRayTelCalorimeterHit.hh"

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class GammaRayTelDetectorConstruction;
class G4HCofThisEvent;
class G4Step;

class GammaRayTelCalorimeterSD: public G4VSensitiveDetector {
public:
	explicit GammaRayTelCalorimeterSD(G4String name);

	~GammaRayTelCalorimeterSD() override;

	void Initialize(G4HCofThisEvent* event) override;

	auto ProcessHits(G4Step *step, G4TouchableHistory *history) -> G4bool  override;

	void EndOfEvent(G4HCofThisEvent *collection) override;

	void clear() override;

	void DrawAll() override;

	void PrintAll() override;

private:
	GammaRayTelCalorimeterHitsCollection *calorimeterCollection;

	GammaRayTelDetectorConstruction *detector;

	G4int *calHitXID;

	G4int *calHitYID;

	G4int numberOfCALLayers;

	G4int numberOfCALBars;

	G4int numberOfCALChannels;
};
#endif
