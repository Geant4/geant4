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
// $Id: HadrontherapyCalorimeterSD.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#ifndef HadrontherapyCalorimeterSD_h
#define HadrontherapyCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "HadrontherapyHit.hh"
class HadrontherapyDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
class HadrontherapyRunAction;
class HadrontherapyHit;

// -------------------------------------------------------------
class HadrontherapyCalorimeterSD : public G4VSensitiveDetector
{
public:
 
  HadrontherapyCalorimeterSD(G4String, HadrontherapyDetectorConstruction*);
  ~HadrontherapyCalorimeterSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void PrintAll();

private:
  G4String filename;
  G4int sliceID[50000];
  G4double energy[50000];
  G4double depth;
  
  HadrontherapyDetectorConstruction* Detector;

  G4double backEnergy;
  G4double leakEnergy;
  G4double delta;
  G4double depthMax;
  G4double tkinold;
  G4bool   part_is_out;
  G4int evno;
  G4int evnOld;
  G4int trIDold;
  G4int NbOfLayer;
  HadrontherapyRunAction* p_Run;
  HadrontherapyHitsCollection* CalCollection;  
};
#endif
























