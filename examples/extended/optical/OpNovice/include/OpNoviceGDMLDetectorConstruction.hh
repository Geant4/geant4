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
#ifndef OpNoviceGDMLDetectorConstruction_h
#define OpNoviceGDMLDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4GDMLParser;
class OpNoviceDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceGDMLDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  OpNoviceGDMLDetectorConstruction(G4String fname);
  ~OpNoviceGDMLDetectorConstruction() override;
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  void ReadGDML();
  void UpdateGeometry();
  void SetDumpGdml(G4bool);
  G4bool IsDumpGdml() const;
  void SetVerbose(G4bool fverbose);
  G4bool IsVerbose() const;
  void SetDumpGdmlFile(G4String fDumpGdmlFile);
  G4String GetDumpGdmlFileName() const;

 private:
  OpNoviceGDMLDetectorConstruction& operator=(
    const OpNoviceGDMLDetectorConstruction& right);
  OpNoviceGDMLDetectorConstruction(const OpNoviceGDMLDetectorConstruction&);

  OpNoviceDetectorMessenger* fDetectorMessenger = nullptr;
  G4GDMLParser* fParser = nullptr;

  G4String fGdmlFile;
  G4String fDumpGdmlFileName = "OpNovice_dump.gdml";
  G4bool fVerbose = false;
  G4bool fDumpGdml = false;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
#endif
