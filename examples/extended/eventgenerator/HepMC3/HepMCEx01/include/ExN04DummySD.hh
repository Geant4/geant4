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
/// \file eventgenerator/HepMC/HepMCEx01/include/ExN04DummySD.hh
/// \brief Definition of the ExN04DummySD class
//
//

// Dummy sensitive used only to flag sensitivity
// in cells of RO geometry.
//

#ifndef ExN04DummySD_h
#define ExN04DummySD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class ExN04DummySD : public G4VSensitiveDetector {
public:
  ExN04DummySD();
  ~ExN04DummySD() {}

  virtual void Initialize(G4HCofThisEvent*) {}
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*) {return false;}
  virtual void EndOfEvent(G4HCofThisEvent*) {}
  virtual void clear() {}
  virtual void DrawAll() {}
  virtual void PrintAll() {}
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04DummySD::ExN04DummySD()
  : G4VSensitiveDetector("dummySD")
{
}

#endif
