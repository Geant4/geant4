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
// This example is provided by the Geant4-DNA collaboration
// dnadamage2 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "ScoreSpecies.hh"
#include "ScoreLET.hh"
#include "ScoreStrandBreaks.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  fDetDir = new G4UIdirectory("/det/");
  fDetDir->SetGuidance("Detector control.");

  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/det/setSize",this);
  fSizeCmd->SetDefaultUnit("um");

  fpOffSetFileUI = new G4UIcmdWithAString("/det/OffSetFile", this);
  fpPlasmidNbUI  = new G4UIcmdWithAnInteger("/det/NbOfPlasmids", this);

  fpPlasmidFile = new G4UIcmdWithAString("/det/PlasmidFile", this);
  fpUseDNA      = new G4UIcmdWithABool("/det/UseDNAVolumes", this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction()
{
  delete fSizeCmd;
  delete fpOffSetFileUI;
  delete fpPlasmidFile;
  delete fpUseDNA;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fSizeCmd)
  {
    G4double size = fSizeCmd->GetNewDoubleValue(newValue);
    SetSize(size);
  }

  if (command == fpPlasmidNbUI)
  {
    fNbOfPlasmids = fpPlasmidNbUI->GetNewIntValue(newValue);
  }

  if (command == fpOffSetFileUI) 
  {
    G4String File = newValue;
    ReadOffsetFile(File);
  }

  if (command == fpPlasmidFile)
  {
    fPlasmidFile = newValue;
  }

  if (command == fpUseDNA)
  {
    fUseDNAVolumes = fpUseDNA->GetNewBoolValue(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetSize(G4double size)
{

  G4GeometryManager::GetInstance()->OpenGeometry();
  fPlasmidEnvelope->SetRadius(size/2);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "#### the geometry has been modified " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Clear DNA information Containers
  fSampleDNANames.clear();
  fSampleDNAPositions.clear();
  fSampleDNADetails.clear();
  fDNANames.clear();
  fDNAPositions.clear();
  fDNADetails.clear();

  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();
  G4Material*   normalWater = man->FindOrBuildMaterial("G4_WATER");
  G4Material*   waterWorld  = man->BuildMaterialWithNewDensity("G4_WATER_WORLD",
                                                               "G4_WATER",
                                                               1.0*g/cm/cm/cm);

  //
  // World
  G4Box* solidWenvelope        = new G4Box("SWorld", 1 * um, 1 * um, 1 * um);
  G4LogicalVolume* logicWorld  = new G4LogicalVolume(solidWenvelope, 
                                                     waterWorld, 
                                                     "LWorld");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, 
                                                   G4ThreeVector(), 
                                                   "PWorld", 
                                                   logicWorld, 
                                                   0, 
                                                   false, 
                                                   0);

  // Molecule
  fPlasmidEnvelope = new G4Orb("plasmidEnvelope", 0.5*fWorldSize);

  G4LogicalVolume* tlogicPlasmid = new G4LogicalVolume(fPlasmidEnvelope, 
                                                       normalWater,
                                                       "PlasmidVolume");
  new G4PVPlacement(0,
                    G4ThreeVector(), 
                    tlogicPlasmid, 
                    "plasmid", 
                    logicWorld, 
                    false, 
                    0);

  // World Vis Attributes
  G4VisAttributes worldVis(G4Colour(0.0, 0.0, 1.0) );
  logicWorld->SetVisAttributes(worldVis);

  PhysGeoImport GeoImport = PhysGeoImport();

  G4LogicalVolume* logicStraightVoxel;

  logicStraightVoxel = GeoImport.CreateLogicVolumeXYZ(fPlasmidFile);

  fSampleDNANames     = GeoImport.GetMoleculesNames();
  fSampleDNAPositions = GeoImport.GetMoleculesPositions();
  fSampleDNADetails   = GeoImport.GetMoleculesDetails();

  for ( int i = 0; i < fNbOfPlasmids; i++ ) {
    if (fUseDNAVolumes)
      new G4PVPlacement(0, 
                        fVOffset[i], 
                        logicStraightVoxel, 
                        "VoxelStraight", 
                        logicWorld, 
                        true, 
                        i);
    AddDNAInformation(i,fVOffset[i]);
  }

  //always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare World as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* mfDetector =
  new G4MultiFunctionalDetector("mfDetector");

  // LET scorer:
  //  - scores restricted or unrestricted LET

  ScoreLET* LET = new ScoreLET("LET");
  mfDetector->RegisterPrimitive(LET);

  // Species scorer:
  //  - scores number of species over time
  //  - compute the radiochemical yields (G values)

  G4VPrimitiveScorer* primitivSpecies = new ScoreSpecies("Species");
  mfDetector->RegisterPrimitive(primitivSpecies);

  // SB Scorer: Requires access to Geometry
  //  - scores number of Direct SB
  //  - scores number of Indirect SB

  G4VPrimitiveScorer* primitiveSB = new ScoreStrandBreaks("StrandBreaks",
                                                         this,&fWorldSize);
  mfDetector->RegisterPrimitive(primitiveSB);

  // Attach Detectors to Plasmid Volumes
  G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for ( size_t t = 0; t < theLogicalVolumes->size(); t++ ) {
    G4String lVolumeName = (*theLogicalVolumes)[t]->GetName();
    if (lVolumeName != "LWorld") {
      (*theLogicalVolumes)[t]->SetSensitiveDetector(mfDetector);
    }
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(mfDetector);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ReadOffsetFile(G4String file) {
  if (fVOffset.size() > 0)
    fVOffset.clear();

  std::ifstream OffSetFile;
  OffSetFile.open(file);
  G4double x, y, z;

  if (!OffSetFile) {
    G4cout << "Plasmid Offset positions file not found!!!" << G4endl;
    exit(1);
  }

  else {
    while (OffSetFile >> x >> y >> z) {
      fVOffset.push_back(G4ThreeVector(x*nm,y*nm,z*nm));
    }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::AddDNAInformation(G4int copy, G4ThreeVector offset) {
  for(size_t i = 0; i < fSampleDNANames.size(); i++) {
    fDNANames.push_back(fSampleDNANames[i]);
    fDNAPositions.push_back(fSampleDNAPositions[i] + offset);

    std::vector<G4int> Details = fSampleDNADetails[i];
    Details[0] = copy;
    fDNADetails.push_back(Details);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
