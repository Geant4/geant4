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
#ifndef TILE05SENSITIVEDETECTOR_H
#define TILE05SENSITIVEDETECTOR_H

using namespace std;
#include "G4VSensitiveDetector.hh"
#include "G4VGFlashSensitiveDetector.hh"
#include "G4GFlashSpot.hh"
#include "ExGflashDetectorConstruction.hh"
#include "ExGflashHitsCollection.hh"
#include "globals.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class ExGflashSensitiveDetector: public G4VSensitiveDetector ,public G4VGFlashSensitiveDetector{
public:
	ExGflashSensitiveDetector(G4String, ExGflashDetectorConstruction* det);
	~ExGflashSensitiveDetector();
	
	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	G4bool ProcessHits(G4GFlashSpot*aSpot,G4TouchableHistory*); 
	void EndOfEvent(G4HCofThisEvent*);
private:
	ExGflashHitsCollection* caloHitsCollection;
	G4double edep;
	ExGflashDetectorConstruction* Detector;
};
#endif
