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
#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* detConstruction)
  : fpDetectorConstruction(detConstruction)
  , fpWorldGeometryDirectory(0)
  , fpWorldSideLength(0)
  , fpCellGeometryDirectory(0)
  , fpCellRadius(0)
{
  // World geometry
  fpWorldGeometryDirectory = new G4UIdirectory("/world/");
  fpWorldGeometryDirectory->SetGuidance("Commands for world geometry params");

  fpWorldSideLength = new G4UIcmdWithADoubleAndUnit("/world/worldSize", this);
  fpWorldSideLength->SetGuidance("Side length for the world");
  fpWorldSideLength->SetParameterName("Side length", false);

  // Cell geometry
  fpCellGeometryDirectory = new G4UIdirectory("/cell/");
  fpCellGeometryDirectory->SetGuidance("Commands for world geometry params");
  fpCellRadius = new G4UIcmdWith3VectorAndUnit("/cell/radiusSize", this);
  fpCellRadius->SetGuidance(
    "Set Semi-Major axes for cell (x, y, z) - unset, whole world is water");
  fpCellRadius->SetParameterName("xradius", "yradius", "zradius", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  // World
  delete fpCellRadius;
  delete fpWorldSideLength;
  delete fpCellGeometryDirectory;
  delete fpWorldGeometryDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  // Geometry
  if(command == fpWorldSideLength)
  {
    fpDetectorConstruction->SetWorldSideLength(
      ((G4UIcmdWithADoubleAndUnit*) command)->GetNewDoubleValue(newValue));
  }
  else if(command == fpCellRadius)
  {
    fpDetectorConstruction->SetCellRadius(
      ((G4UIcmdWith3VectorAndUnit*) command)->GetNew3VectorValue(newValue));
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
