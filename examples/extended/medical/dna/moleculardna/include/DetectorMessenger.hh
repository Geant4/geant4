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

#ifndef MOLECULAR_DETECTOR_MESSENGER_HH
#define MOLECULAR_DETECTOR_MESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;

class G4UIcmdWith3VectorAndUnit;

class G4UIcmdWithADoubleAndUnit;

class G4UIdirectory;

class G4UIcommand;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger : public G4UImessenger
{
 public:
  explicit DetectorMessenger(DetectorConstruction*);

  ~DetectorMessenger() override;

  void SetNewValue(G4UIcommand*, G4String) override;

 private:
  DetectorConstruction* fpDetectorConstruction;

  // Related to geometry
  G4UIdirectory* fpWorldGeometryDirectory;
  G4UIcmdWithADoubleAndUnit* fpWorldSideLength;

  G4UIdirectory* fpCellGeometryDirectory;
  G4UIcmdWith3VectorAndUnit* fpCellRadius;
};

#endif  // MOLECULAR_DETECTOR_MESSENGER_HH
