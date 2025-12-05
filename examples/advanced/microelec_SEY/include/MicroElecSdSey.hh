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
//---------------------------------------------------------------------------
//
// ClassName:   MicroElecSdSey
//
// Description: The process to kill e- to save CPU
//
// Author:      C. Inguimbert 16/02/2022
//				ONERA
//----------------------------------------------------------------------------
//
// Class description:
//
// SEY : Secondary Electron Emission Yield
// detecteor to be used to count the number of secondary electrons
// emitted by a surface of irradiated material
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#ifndef SD_MicroElecSdSey
#define SD_MicroElecSdSey 1

#include	<iostream>
#include	<fstream>
#include	<vector>

#include	"MicroElecHitSey.hh"
#include	"globals.hh"
#include	"G4Timer.hh"
#include	"G4VSensitiveDetector.hh"
#include	"G4GeneralParticleSource.hh"
#include	"G4Vector3D.hh"
#include	"G4LogicalVolume.hh"


class G4Step;
class G4HCofThisEvent;

class MicroElecSdSey : public G4VSensitiveDetector
{
	
public:

	MicroElecSdSey(const G4String& name,
		const G4String& hitsCollectionName);
	~MicroElecSdSey();

	G4double nbHitsTotal_;
	G4double nbHitsTotal_bis;

	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	void EndOfEvent(G4HCofThisEvent*);

	G4double GetCompteurPrim() { return compteurPrimaire; };
	G4double GetCompteurSec() { return compteurSec; };
	G4double GetCompteurTot() { return compteurTot; };
	G4double GetCompteur50() { return compteur50; };

	void ResetCounters() {
		compteurPrimaire = 0; compteurSec = 0; compteurTot = 0; compteur50=0;
	}

	inline void SetEdep(G4double de) { Sensitive_Detector_Edep = de; };
	inline G4double GetEdep() { return Sensitive_Detector_Edep; };
	inline void AddEdep(G4double de) { Sensitive_Detector_Edep = Sensitive_Detector_Edep + de; };

private:
	G4double compteurPrimaire, compteurSec, compteurTot, compteur50, nbPrim, nbSec, nbSup50;
	MicroElecHitSeyCollection* fHitsCollection;
	

public:
	G4bool FileFlag;
	G4double Sensitive_Detector_Edep;

};

#endif
