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
/// file:DNAWorld.hh
/// brief :
/*
 * This class contains the physical placements of DNA volumes.
 * In order to keep all the volume placement in the same place, volumes
 * are defined and placed in DNAGeometry, and their parent volume
 * is passed to this class via SetDNAVolumePointer()
 *
 * This volume is then placed in the world volume.
 * It is important that a physics process is created to ensure that this
 * world functions as a layered mass geometry.
 *
 * The aim of this function is to allow the chemistry to run without seeing
 * DNA volumes, in the geometry, but still allowing the chemistry to access
 * the volume positions through an alternate data structure.
 */

#ifndef MOLECULAR_DNA_WORLD_HH
#define MOLECULAR_DNA_WORLD_HH

#include "G4VUserParallelWorld.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAWorld : public G4VUserParallelWorld
{
 public:
  DNAWorld();

  ~DNAWorld() override;

  void Construct() override;

  void SetDNAVolumePointer(G4LogicalVolume* lv) { fpDNAVolumePointer = lv; };

  void SetDNAVolumeTranslation(G4ThreeVector* t)
  {
    fpDNAVolumeTranslation = t;
  };

  void SetDNAVolumeRotation(G4RotationMatrix* rot)
  {
    fpDNAVolumeRotation = rot;
  };

 private:
  G4LogicalVolume* fpDNAVolumePointer = nullptr;
  G4ThreeVector* fpDNAVolumeTranslation = nullptr;
  G4RotationMatrix* fpDNAVolumeRotation = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_DNA_WORLD_HH
