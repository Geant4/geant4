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
//      ------------ GammaRayTelAnalysis  ------
//           by R.Giannitrapani, F. Longo & G.Santin (30 nov 2000)
//
// 03.04.2013 F.Longo and L.Pandola
// - migrated to G4AnalysisManager
//
// 07.12.2001 A.Pfeiffer
// - integrated Guy's addition of the ntuple
//
// 06.12.2001 A.Pfeiffer
// - updating to new design (singleton)
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
// -------------------------------------------------------------------
// Class description:
// Example of analysis in a simulation application (histograms, ntuples etc.)
// This class follows the singleton design pattern; 
// it is responsible for the analysis management and algorithms 
//
// -------------------------------------------------------------------//

#ifndef GammaRayTelAnalysis_h 
#define GammaRayTelAnalysis_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "G4AnalysisManager.hh"

class GammaRayTelAnalysisMessenger;
class GammaRayTelDetectorConstruction;

class GammaRayTelAnalysis {
public:
	~GammaRayTelAnalysis();

	void BeginOfRun();

	void EndOfRun();

	void EndOfEvent(G4int flag);

	void Init();

	void Finish();

	void SetHisto2DMode(G4String value) {
		histo2DMode = value;
	}

	auto GetHisto2DMode() -> G4String {
		return histo2DMode;
	}

	void InsertPositionXZ(G4double x, G4double z);

	void InsertPositionYZ(G4double y, G4double z);

	void InsertEnergy(G4double energy);

	void InsertHits(G4int planeNumber);

	void setNtuple(G4double energy, G4int planeNumber, G4double x, G4double y, G4double z);

	static auto getInstance() -> GammaRayTelAnalysis*;

private:
	GammaRayTelAnalysis();

	void Plot();

	static GammaRayTelAnalysis *instance;

	const GammaRayTelDetectorConstruction *detector{nullptr};

	G4String histo2DMode;

	G4String histogramFileName;

	GammaRayTelAnalysisMessenger *analysisMessenger;
};
#endif
