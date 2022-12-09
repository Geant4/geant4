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
//      ------------ GammaRayTelTrackerSD  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelTrackerSD_h
#define GammaRayTelTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class GammaRayTelDetectorConstruction;
class G4HCofThisEvent;
class G4Step;

#include "GammaRayTelTrackerHit.hh"

class GammaRayTelTrackerSD: public G4VSensitiveDetector {
public:
    explicit GammaRayTelTrackerSD(G4String name);

	~GammaRayTelTrackerSD() override;

	void Initialize(G4HCofThisEvent *event) override;

	auto ProcessHits(G4Step *step, G4TouchableHistory *history) -> G4bool override;

	void EndOfEvent(G4HCofThisEvent *collection) override;

	void clear() override;

	void DrawAll() override;

	void PrintAll() override;

private:
	GammaRayTelTrackerHitsCollection *trackerCollection;

	GammaRayTelDetectorConstruction *detector;

	G4int *tkrHitXID;

	G4int *tkrHitYID;

	G4int numberOfTKRLayers;

	G4int numberOfTKRStrips;

	G4int numberOfTKRChannels;
};
#endif
