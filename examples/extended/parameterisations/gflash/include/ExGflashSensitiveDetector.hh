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
// $Id: ExGflashSensitiveDetector.hh 94396 2015-11-13 13:37:16Z gcosmo $
//
/// \file parameterisations/gflash/include/ExGflashSensitiveDetector.hh
/// \brief Definition of the ExGflashSensitiveDetector class
//
#ifndef TILE05SENSITIVEDETECTOR_H
#define TILE05SENSITIVEDETECTOR_H

#include "G4VSensitiveDetector.hh"
#include "G4VGFlashSensitiveDetector.hh"
#include "ExGflashHitsCollection.hh"
#include "globals.hh"

class G4GFlashSpot;
class ExGflashDetectorConstruction;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ExGflashSensitiveDetector: public G4VSensitiveDetector, 
                                 public G4VGFlashSensitiveDetector {
public:
  ExGflashSensitiveDetector(G4String, ExGflashDetectorConstruction* det);
  ~ExGflashSensitiveDetector();
  
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  virtual G4bool ProcessHits(G4GFlashSpot*aSpot,G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);

private:
  ExGflashHitsCollection* fCaloHitsCollection;
  ExGflashDetectorConstruction* fDetector;
  G4int fHCID;
};

#endif
