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
/// file:DNAWorld.cc
/// brief:

#include "DNAWorld.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAWorld::DNAWorld()
  : G4VUserParallelWorld("DNAWorld")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAWorld::~DNAWorld() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAWorld::Construct()
{
  G4VPhysicalVolume* ghostWorld = GetWorld();
  if(fpDNAVolumePointer == nullptr)
  {
    G4cout << "Error: The DNA volume pointer is not set in the DNAWorld"
           << G4endl;
  }
  else
  {
    if(fpDNAVolumeTranslation == nullptr)
    {
      fpDNAVolumeTranslation = new G4ThreeVector(0., 0., 0.);
    }
    new G4PVPlacement(fpDNAVolumeRotation, *fpDNAVolumeTranslation, "DNA",
                      fpDNAVolumePointer, ghostWorld, false, 0, false);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
