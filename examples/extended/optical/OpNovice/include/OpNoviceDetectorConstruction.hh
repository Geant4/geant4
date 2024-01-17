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
/// \file OpNovice/include/OpNoviceDetectorConstruction.hh
/// \brief Definition of the OpNoviceDetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef OpNoviceDetectorConstruction_h
#define OpNoviceDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <CLHEP/Units/SystemOfUnits.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceDetectorMessenger;

class OpNoviceDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  OpNoviceDetectorConstruction();
  ~OpNoviceDetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  void SetVerbose(G4bool verbose);
  G4bool IsVerbose() const;

  void SetDumpGdml(G4bool);
  G4bool IsDumpGdml() const;
  void SetDumpGdmlFile(G4String);
  G4String GetDumpGdmlFile() const;

 private:
  void PrintError(G4String);

  OpNoviceDetectorMessenger* fDetectorMessenger = nullptr;
  G4String fDumpGdmlFileName = "OpNovice_dump.gdml";

  G4double fWorld_x = 15.*CLHEP::m;
  G4double fWorld_y = 15.*CLHEP::m;
  G4double fWorld_z = 15.*CLHEP::m;

  G4double fExpHall_x = 10.*CLHEP::m;
  G4double fExpHall_y = 10.*CLHEP::m;
  G4double fExpHall_z = 10.*CLHEP::m;

  G4double fTank_x = 5.*CLHEP::m;
  G4double fTank_y = 5.*CLHEP::m;
  G4double fTank_z = 5.*CLHEP::m;

  G4double fBubble_x = 0.5*CLHEP::m;
  G4double fBubble_y = 0.5*CLHEP::m;
  G4double fBubble_z = 0.5*CLHEP::m;
  
  G4bool fVerbose = false;
  G4bool fDumpGdml = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpNoviceDetectorConstruction_h*/
