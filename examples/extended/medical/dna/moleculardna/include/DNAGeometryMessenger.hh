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
/// \file DNAGeometryMessenger.hh

#ifndef MOLECULAR_DNA_MESSENGER_HH
#define MOLECULAR_DNA_MESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

#include <memory>

class DNAGeometry;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAGeometryMessenger : public G4UImessenger
{
 public:
  explicit DNAGeometryMessenger(DNAGeometry*);

  ~DNAGeometryMessenger() override = default;

  void SetNewValue(G4UIcommand*, G4String) override;

 protected:
 private:
  DNAGeometry* fpDNAGeometry;

  // Related to "voxel" placements
  std::unique_ptr<G4UIdirectory> fpDNADirectory;
  std::unique_ptr<G4UIcmdWithAString>
    fpNewPlacementVolumeFile;  // specifies molecule locations
  std::unique_ptr<G4UIcmdWithAString>
    fpVoxelPlacementsFile;  // specifies voxel locations
  std::unique_ptr<G4UIcmdWith3VectorAndUnit> fpVoxelSideLength;  // voxel size
  std::unique_ptr<G4UIcmdWith3VectorAndUnit>
    fpFractalScaling;  // multiplier for de-dimensionalised voxel locations
  std::unique_ptr<G4UIcmdWithABool>
    fpAnglesAsPi;  // take voxel placement angles as multiples of pi
  std::unique_ptr<G4UIcmdWithABool>
    fpCustomMoleculeSizes;  // turn off default molecule sizes and use custom
  std::unique_ptr<G4UIcmdWithAString> fpAddMoleculeSize;  // set a molecule size

  // Related to processing
  std::unique_ptr<G4UIcmdWithABool> fpCheckDNAOverlaps;
  std::unique_ptr<G4UIcmdWithAnInteger> fpVerbosity;
  std::unique_ptr<G4UIcmdWithAnInteger> fpSmartless;

  // damage-like
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpIntRangeDirect;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpRadicalKillDistance;
  std::unique_ptr<G4UIcmdWithABool> fpUseHistoneScav;

  // vis
  std::unique_ptr<G4UIcmdWithABool> fpDrawCellVolumes;
  // testing
  std::unique_ptr<G4UIdirectory> fpTestDirectory;
  std::unique_ptr<G4UIcmdWithoutParameter> fpChromosomeTest;
  std::unique_ptr<G4UIcmdWithoutParameter> fpBasePairTest;
  std::unique_ptr<G4UIcmdWithoutParameter> fpUniqueIDTest;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_DNA_MESSENGER_HH
