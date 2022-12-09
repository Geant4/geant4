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
#include "DNAGeometryMessenger.hh"
#include "DNAGeometry.hh"
#include "UtilityFunctions.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAGeometryMessenger::DNAGeometryMessenger(DNAGeometry* dnaGeometry)
  : fpDNAGeometry(dnaGeometry)
  , fpDNADirectory(new G4UIdirectory("/dnageom/"))
  , fpNewPlacementVolumeFile(
      new G4UIcmdWithAString("/dnageom/placementVolume", this))
  , fpVoxelPlacementsFile(
      new G4UIcmdWithAString("/dnageom/definitionFile", this))
  , fpVoxelSideLength(
      new G4UIcmdWith3VectorAndUnit("/dnageom/placementSize", this))
  , fpFractalScaling(
      new G4UIcmdWith3VectorAndUnit("/dnageom/fractalScaling", this))
  , fpAnglesAsPi(new G4UIcmdWithABool(
      "/dnageom/setVoxelPlacementAnglesAsMultiplesOfPi", this))
  , fpCustomMoleculeSizes(
      new G4UIcmdWithABool("/dnageom/useCustomMoleculeSizes", this))
  , fpAddMoleculeSize(new G4UIcmdWithAString("/dnageom/moleculeSize", this))
  , fpCheckDNAOverlaps(new G4UIcmdWithABool("/dnageom/checkOverlaps", this))
  , fpVerbosity(new G4UIcmdWithAnInteger("/dnageom/verbose", this))
  , fpSmartless(new G4UIcmdWithAnInteger("/dnageom/setSmartVoxels", this))
  , fpIntRangeDirect(
      new G4UIcmdWithADoubleAndUnit("/dnageom/interactionDirectRange", this))
  , fpRadicalKillDistance(
      new G4UIcmdWithADoubleAndUnit("/dnageom/radicalKillDistance", this))
  , fpUseHistoneScav(
      new G4UIcmdWithABool("/dnageom/activateHistoneScavenging", this))
  , fpDrawCellVolumes(new G4UIcmdWithABool("/dnageom/drawCellVolumes", this))
  , fpTestDirectory(new G4UIdirectory("/dnatests/"))
  , fpChromosomeTest(new G4UIcmdWithoutParameter("/dnatests/chromosome", this))
  , fpBasePairTest(new G4UIcmdWithoutParameter("/dnatests/basepairs", this))
  , fpUniqueIDTest(new G4UIcmdWithoutParameter("/dnatests/uniqueid", this))
{
  // Fractals
  fpDNADirectory->SetGuidance("Commands to control the fractal.");

  fpNewPlacementVolumeFile->SetGuidance("Set a placement volume");
  fpNewPlacementVolumeFile->SetGuidance("format: name path");
  fpNewPlacementVolumeFile->SetParameterName("name path", false);

  fpVoxelPlacementsFile->SetGuidance(
    "Path to file that defines placement locations");
  fpVoxelPlacementsFile->SetParameterName("path", false);

  fpVoxelSideLength->SetGuidance("Side length for each placement (x, y, z)");
  fpVoxelSideLength->SetParameterName("xlength", "ylength", "zlength", false);

  fpFractalScaling->SetGuidance("Scaling of XYZ in fractal definition file");
  fpFractalScaling->SetParameterName("xlength", "ylength", "zlength", false);

  fpAnglesAsPi->SetGuidance(
    "Take the angles in the voxel placement file as multiples of pi");
  fpAnglesAsPi->SetGuidance(
    "E.g. set to true if the angle 0.5 should mean 90 degrees");
  fpAnglesAsPi->SetParameterName("true/false angles as pi", false);
  fpAnglesAsPi->SetDefaultValue(false);

  fpAnglesAsPi->SetGuidance("Enable custom molecule sizes");
  fpAnglesAsPi->SetGuidance("These can now be set via /dnageom/moleculeSize");
  fpAnglesAsPi->SetParameterName("true/false custom sizes", false);
  fpAnglesAsPi->SetDefaultValue(false);

  fpAddMoleculeSize->SetGuidance("Set a molecule size in angstrom.");
  fpAddMoleculeSize->SetGuidance("format: molecule_name x y z");
  fpAddMoleculeSize->SetGuidance("E.G.: PHOSPHATE 3 4 5");
  fpAddMoleculeSize->SetGuidance("Note: molecule names are case insensitive");
  fpAddMoleculeSize->SetParameterName("name x y z units", false);

  // control
  fpCheckDNAOverlaps->SetGuidance("Check overlaps in DNA geometry region");
  fpCheckDNAOverlaps->SetParameterName("true/false check overlaps", false);

  fpVerbosity->SetGuidance("Verbosity for DNA geometry");
  fpVerbosity->SetParameterName("int verbose level", false);

  fpSmartless->SetGuidance("Optimisation value (int) for smart voxels");
  fpSmartless->SetGuidance("The G4 default is 2");
  fpSmartless->SetParameterName("Optimasation value", false);

  fpIntRangeDirect->SetGuidance(
    "Critical range to start recording SSBs from direct effects");
  fpIntRangeDirect->SetParameterName("Range", false);

  fpRadicalKillDistance->SetGuidance(
    "Distance from base pairs at which to kill radicals");
  fpRadicalKillDistance->SetParameterName("Range", false);

  fpUseHistoneScav->SetGuidance(
    "Activate Histone scavenging function with default radius");
  fpUseHistoneScav->SetGuidance(
    "Radius can be controlled with /dnageom/histoneScavengingRadius");
  fpUseHistoneScav->SetParameterName("true/false, set histone scavenging on",
                                     false);

  fpDrawCellVolumes->SetGuidance(
    "Draw cell/chromosome volumes rather than DNA (makes DNA invisible)");
  fpDrawCellVolumes->SetParameterName("true/false draw cell volumes", false);

  // Tests
  fpTestDirectory->SetGuidance("Tests of the DNA geometry");

  fpChromosomeTest->SetGuidance("Test Chromosome Placement classes");

  fpBasePairTest->SetGuidance("Test Base Pair indices are correct");

  fpUniqueIDTest->SetGuidance("Test Unique ID Algorithm");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  // Fractals
  if(command == fpNewPlacementVolumeFile.get())
  {
    std::vector<G4String> cmd = utility::Split(newValue, ' ');
    // input validation
    if(!(cmd.size() == 2 || cmd.size() == 3))
    {
      G4cout << "Invalid input. Command takes 2 string parameters only"
             << G4endl << "separated by a space, or three parameters" << G4endl
             << "where the last is a boolean." << G4endl;
    }
    else
    {
      G4String vname = (G4String) cmd[0];
      G4String path  = (G4String) cmd[1];
      G4bool twist   = false;
      if(cmd.size() == 3)
      {
        if((G4String) cmd[2] == "true")
        {
          twist = true;
        }
        else if((G4String) cmd[2] != "false")
        {
          G4ExceptionDescription errmsg;
          errmsg << "Invalid input. Third value must be a boolean"
                 << "written as true or false." << G4endl;
          G4Exception("DNAGeometryMessenger::SetNewValue", "DNAGeometry001",
                      FatalException, errmsg);
        }
      }

      if(utility::Path_exists(path))
      {
        fpDNAGeometry->AddVoxelFile(vname, path, twist);
      }
      else
      {
        G4ExceptionDescription errmsg;
        errmsg << path << " has wrong path or does not exist :";
        G4Exception("DNAGeometryMessenger::SetNewValue", "DNAGeometry002",
                    FatalException, errmsg);
      }
    }
  }
  else if(command == fpVoxelPlacementsFile.get())
  {
    G4String path = newValue;
    if(utility::Path_exists(path))
    {
      fpDNAGeometry->SetFractalFilename(path);
    }
    else
    {
      G4ExceptionDescription errmsg;
      errmsg << "The fractal file path has wrong path or does not exist : "
             <<path<< G4endl;
      G4Exception("DNAGeometryMessenger::SetNewValue", "DNAGeometry003",
                  FatalException, errmsg);
    }
  }
  else if(command == fpVoxelSideLength.get())
  {
    fpDNAGeometry->SetVoxelSideLength(
      G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
  }
  else if(command == fpFractalScaling.get())
  {
    fpDNAGeometry->SetFractalScaling(
      G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(newValue));
  }
  else if(command == fpAnglesAsPi.get())
  {
    fpDNAGeometry->SetFractalAnglesAsPi(
      G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
  else if(command == fpCustomMoleculeSizes.get())
  {
    fpDNAGeometry->EnableCustomMoleculeSizes(
      G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
  else if(command == fpAddMoleculeSize.get())
  {
    std::vector<G4String> cmd = utility::Split(newValue, ' ');
    // input validation
    if(cmd.size() == 4)
    {
      G4cout << "Invalid input. Command takes 4  parameters only" << G4endl
             << "separated by a space, e.g. phosphate 1 2 3" << G4endl;
    }
    else
    {
      try
      {
        G4String name               = cmd[0];
        G4int x_size                = std::stod(cmd[1]);
        G4int y_size                = std::stod(cmd[2]);
        G4int z_size                = std::stod(cmd[3]);
        G4ThreeVector molecule_size = G4ThreeVector(x_size, y_size, z_size);
        fpDNAGeometry->AddChangeMoleculeSize(name, molecule_size);
      } catch(const std::invalid_argument& ia)
      {
        G4cerr << "Invalid argument to convert to double: " << ia.what()
               << G4endl;
      }
    }
  }
  // control
  else if(command == fpCheckDNAOverlaps.get())
  {
    fpDNAGeometry->SetOverlaps(G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
  else if(command == fpVerbosity.get())
  {
    fpDNAGeometry->SetVerbosity(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
  else if(command == fpSmartless.get())
  {
    fpDNAGeometry->SetSmartless(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
  else if(command == fpIntRangeDirect.get())
  {
    fpDNAGeometry->SetDirectInteractionRange(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
  }
  else if(command == fpRadicalKillDistance.get())
  {
    fpDNAGeometry->SetRadicalKillDistance(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
  }
  else if(command == fpUseHistoneScav.get())
  {
    fpDNAGeometry->SetHistoneScav(G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
  else if(command == fpDrawCellVolumes.get())
  {
    fpDNAGeometry->SetDrawCellVolumes(
      G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
  else if(command == fpChromosomeTest.get())
  {
    fpDNAGeometry->ChromosomeTest();
  }
  else if(command == fpBasePairTest.get())
  {
    fpDNAGeometry->BasePairIndexTest();
  }
  else if(command == fpUniqueIDTest.get())
  {
    fpDNAGeometry->UniqueIDTest();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
