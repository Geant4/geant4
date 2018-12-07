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
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fpMessenger(0),
  fCheckOverlaps(false)
{
  //Select a PDB file name by default
  //otherwise modified by the LoadPDBfile messenger
  //
  fPdbFileName=G4String("1ZBB.pdb");
  fPdbFileStatus=0;
  fChosenOption=11;

  fpDefaultMaterial=0;
  fpWaterMaterial=0;
  fpMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  fChosenOption=11;//Draw atomic with bounding box (visualisation by default)
  fPdbFileStatus=0;//There is no PDB file loaded at this stage

  //Define materials and geometry
  G4VPhysicalVolume* worldPV;
  worldPV=DefineVolumes(fPdbFileName,fChosenOption);
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
  //[G4_WATER]
  //
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_WATER", fromIsotopes);
  fpWaterMaterial = G4Material::GetMaterial("G4_WATER");

  //[Vacuum]
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                 kStateGas, 2.73*kelvin, 3.e-18*pascal);
  fpDefaultMaterial = G4Material::GetMaterial("Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CheckMaterials()
{
  if(!fpDefaultMaterial)
    G4Exception("DetectorConstruction::CheckMaterials",
                "DEFAULT_MATERIAL_NOT_INIT_1",
                FatalException,
                "Default material not initialized.");
                
  if(!fpWaterMaterial)
    G4Exception("DetectorConstruction::CheckMaterials",
                "WATER_MATERIAL_NOT_INIT",
                FatalException,
                "Water material not initialized.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{
  // Geometry parameters
  G4double worldSize  = 1000*1*angstrom;

  if ( !fpDefaultMaterial )
  {
     G4Exception("DetectorConstruction::ConstructWorld",
                "DEFAULT_MATERIAL_NOT_INIT_2",
                FatalException,
                "Default material not initialized.");
  }

  //
  // World
  //
  G4VSolid* worldS
  = new G4Box("World",           // its name
              worldSize/2, worldSize/2, worldSize/2); // its size

  G4LogicalVolume*
  worldLV
  = new G4LogicalVolume(
      worldS,           // its solid
      fpDefaultMaterial,  // its material
      "World");         // its name

  G4VisAttributes *MyVisAtt_ZZ = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));
  MyVisAtt_ZZ ->SetVisibility (false);
  worldLV->SetVisAttributes(MyVisAtt_ZZ);

  G4VPhysicalVolume*
  worldPV
  = new G4PVPlacement(
      0,                // no rotation
      G4ThreeVector(),  // at (0,0,0)
      worldLV,          // its logical volume
      "World",          // its name
      0,                // its mother  volume
      false,            // no boolean operation
      0,                // copy number
      true);            // checking overlaps forced to YES

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes(G4String filename,
    unsigned short int option)
{
  //Clean old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Define Materials
  //
  ConstructMaterials();

  //Reconstruct world not to superimpose on geometries
  G4VPhysicalVolume* worldPV;
  G4LogicalVolume* worldLV;
  worldPV=ConstructWorld();
  worldLV=worldPV->GetLogicalVolume();

  //List of molecules inside pdb file separated by TER keyword
  fpMoleculeList=NULL;
  fpBarycenterList=NULL;

  //'fPDBlib.load' is currently done for each 'DefineVolumes' call.
  int verbosity=0;  //To be implemented
  unsigned short int isDNA;
  fpMoleculeList = fPDBlib.Load(filename,isDNA,verbosity);
  if (fpMoleculeList!=NULL && isDNA==1)
  {
    fPDBlib.ComputeNbNucleotidsPerStrand(fpMoleculeList );
    fpBarycenterList=fPDBlib.ComputeNucleotideBarycenters(fpMoleculeList );
    G4cout<<"This PDB file is DNA."<<G4endl;
  }

  if (fpMoleculeList!=NULL)
    {
    fPdbFileStatus=1;
    }

  if (option==1)
  {
    AtomisticView( worldLV,fpMoleculeList,1.0);
  }
  else
    if (option==2)
    {
      BarycenterView(worldLV,fpBarycenterList);
    }
    else
      if (option==3)
      {
        ResiduesView(worldLV,fpBarycenterList);
      }
      else
        if (option==10)
        {
          DrawBoundingVolume( worldLV,fpMoleculeList);
        }
        else
          if (option==11)
          {
            AtomisticView( worldLV,fpMoleculeList,1.0);
            DrawBoundingVolume( worldLV,fpMoleculeList);
          }
          else
            if (option==12)
            {
              BarycenterView(worldLV,fpBarycenterList);
              DrawBoundingVolume( worldLV,fpMoleculeList);
            }
            else
              if (option==13)
              {
                ResiduesView(worldLV,fpBarycenterList);
                DrawBoundingVolume( worldLV,fpMoleculeList);
              }
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PDBlib DetectorConstruction::GetPDBlib()
{
  return fPDBlib;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Barycenter *DetectorConstruction::GetBarycenterList()
{
  return fpBarycenterList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule *DetectorConstruction::GetMoleculeList()
{
  return fpMoleculeList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////
///////////////// BEGIN atomistic representation
//
void DetectorConstruction::AtomisticView(G4LogicalVolume* worldLV,
    Molecule *moleculeListTemp, double atomSizeFactor)
{
  CheckMaterials();

  // All solids are defined for different color attributes :
  G4double sphereSize  = atomSizeFactor*1*angstrom;
  G4VSolid* atomS_H = new G4Orb("Sphere", sphereSize*1.2);
  G4VSolid* atomS_C = new G4Orb("Sphere", sphereSize*1.7);
  G4VSolid* atomS_O = new G4Orb("Sphere", sphereSize*1.52);
  G4VSolid* atomS_N = new G4Orb("Sphere", sphereSize*1.55);
  G4VSolid* atomS_S = new G4Orb("Sphere", sphereSize*1.8);
  G4VSolid* atomS_P = new G4Orb("Sphere", sphereSize*1.8);
  G4VSolid* atomS_X = new G4Orb("Sphere", sphereSize);  //Undefined

  //Logical volumes and color table for atoms :
  G4LogicalVolume* atomLV_H = new G4LogicalVolume(
      atomS_H, fpWaterMaterial,  "atomLV_H");  // its solid, material, name
  G4VisAttributes * MyVisAtt_H = new G4VisAttributes(
      G4Colour(G4Colour::White()));
  MyVisAtt_H->SetForceSolid(true);
  atomLV_H->SetVisAttributes(MyVisAtt_H);

  G4LogicalVolume* atomLV_C = new G4LogicalVolume(
      atomS_C, fpWaterMaterial, "atomLV_C");  // its solid, material, name
  G4VisAttributes * MyVisAtt_C = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));//Black in CPK, but bad rendering
  MyVisAtt_C->SetForceSolid(true);
  atomLV_C->SetVisAttributes(MyVisAtt_C);

  G4LogicalVolume* atomLV_O = new G4LogicalVolume(
      atomS_O, fpWaterMaterial, "atomLV_O"); // its solid, material, name
  G4VisAttributes * MyVisAtt_O = new G4VisAttributes(
      G4Colour(G4Colour::Red()));
  MyVisAtt_O->SetForceSolid(true);
  atomLV_O->SetVisAttributes(MyVisAtt_O);

  G4LogicalVolume* atomLV_N = new G4LogicalVolume(
      atomS_N, fpWaterMaterial, "atomLV_N"); // its solid, material, name
  G4VisAttributes * MyVisAtt_N = new G4VisAttributes(
      G4Colour(G4Colour(0.,0.,0.5)));//Dark blue
  MyVisAtt_N->SetForceSolid(true);
  atomLV_N->SetVisAttributes(MyVisAtt_N);

  G4LogicalVolume* atomLV_S = new G4LogicalVolume(
      atomS_S, fpWaterMaterial, "atomLV_S"); // its solid, material, name
  G4VisAttributes * MyVisAtt_S = new G4VisAttributes(G4Colour(
      G4Colour::Yellow()));
  MyVisAtt_S->SetForceSolid(true);
  atomLV_S->SetVisAttributes(MyVisAtt_S);

  G4LogicalVolume* atomLV_P = new G4LogicalVolume(
      atomS_P, fpWaterMaterial, "atomLV_P"); // its solid, material, name
  G4VisAttributes * MyVisAtt_P = new G4VisAttributes(
      G4Colour(G4Colour(1.0,0.5,0.)));//Orange
  MyVisAtt_P->SetForceSolid(true);
  atomLV_P->SetVisAttributes(MyVisAtt_P);

  G4LogicalVolume* atomLV_X = new G4LogicalVolume(
      atomS_X, fpWaterMaterial, "atomLV_X"); // its solid, material, name
  G4VisAttributes * MyVisAtt_X = new G4VisAttributes(
      G4Colour(G4Colour(1.0,0.75,0.8)));//Pink, other elements in CKP
  MyVisAtt_X->SetForceSolid(true);
  atomLV_X->SetVisAttributes(MyVisAtt_X);

  //Placement and physical volume construction from memory
  Residue *residueListTemp;
  Atom *AtomTemp;

  int nbAtomTot=0;  //Number total of atoms
  int nbAtomH=0, nbAtomC=0, nbAtomO=0, nbAtomN=0, nbAtomS=0, nbAtomP=0;
  int nbAtomX=0;
  int k=0;

  while (moleculeListTemp)
  {
    residueListTemp = moleculeListTemp->GetFirst();

    k++;
    int j=0;
    while (residueListTemp)
    {
      AtomTemp=residueListTemp->GetFirst();
      j++;

      int startFrom=0;
      int upTo=residueListTemp->fNbAtom; //case Base+sugar+phosphat (all atoms)
      for (int i=0 ; i < (upTo + startFrom) ; i++) //Phosphat or Sugar or Base
      {

        if (AtomTemp->fElement.compare("H") == 0)
        {
          nbAtomH++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_H,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("C") == 0)
        {
          nbAtomC++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_C,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("O") == 0)
        {
          nbAtomO++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_O,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("N") == 0)
        {
          nbAtomN++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_N,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("S") == 0)
        {
          nbAtomS++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_S,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("P") == 0)
        {
          nbAtomP++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_P,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else
        {
          nbAtomX++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_X,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }

        nbAtomTot++;
        //}//End if (i>=startFrom)
        AtomTemp=AtomTemp->GetNext();
      }//end of for (  i=0 ; i < residueListTemp->nbAtom ; i++)

      residueListTemp=residueListTemp->GetNext();
    }//end of while (residueListTemp)

    moleculeListTemp=moleculeListTemp->GetNext();
  }//end of while (moleculeListTemp)

  G4cout << "**************** atomisticView(...) ****************" <<G4endl;
  G4cout << "Number of loaded chains = " << k <<G4endl;
  G4cout << "Number of Atoms      = " << nbAtomTot <<G4endl;
  G4cout << "Number of Hydrogens  = " << nbAtomH <<G4endl;
  G4cout << "Number of Carbons    = " << nbAtomC <<G4endl;
  G4cout << "Number of Oxygens    = " << nbAtomO <<G4endl;
  G4cout << "Number of Nitrogens  = " << nbAtomN <<G4endl;
  G4cout << "Number of Sulfurs    = " << nbAtomS <<G4endl;
  G4cout << "Number of Phosphorus = " << nbAtomP <<G4endl;
  G4cout << "Number of undifined atoms =" << nbAtomX <<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////
///////////////// BEGIN representation for nucleotide barycenters
//
void DetectorConstruction::BarycenterView(G4LogicalVolume* worldLV,
    Barycenter *barycenterListTemp)
{
  CheckMaterials();

  G4VSolid* atomS_ZZ;
  G4LogicalVolume* atomLV_ZZ;
  G4VisAttributes* MyVisAtt_ZZ;

  int k=0;

  while (barycenterListTemp)
  {
    k++;

    atomS_ZZ = new G4Orb("Sphere", (barycenterListTemp->GetRadius())*angstrom);
    atomLV_ZZ = new G4LogicalVolume(atomS_ZZ,fpWaterMaterial,"atomLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Magenta()));
    MyVisAtt_ZZ->SetForceSolid(true);
    atomLV_ZZ->SetVisAttributes(MyVisAtt_ZZ);

    new G4PVPlacement(0, G4ThreeVector(
        barycenterListTemp->fCenterX*1*angstrom,
        barycenterListTemp->fCenterY*1*angstrom,
        barycenterListTemp->fCenterZ*1*angstrom),
        atomLV_ZZ,
        "atomZZ",
        worldLV,
        false,
        0,
        fCheckOverlaps);

    barycenterListTemp=barycenterListTemp->GetNext();
  }//end of while (moleculeListTemp)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////
///////////////// BEGIN representation for Base Sugar Phosphate barycenters
//
void DetectorConstruction::ResiduesView(G4LogicalVolume* worldLV,
    Barycenter *barycenterListTemp)
{
  CheckMaterials();
  G4VisAttributes* MyVisAtt_ZZ;

  G4VSolid* tubS1_ZZ;
  G4LogicalVolume* tubLV1_ZZ;
  G4VSolid* tubS2_ZZ;
  G4LogicalVolume* tubLV2_ZZ;

  G4VSolid* AS_ZZ;
  G4LogicalVolume* ALV_ZZ;
  G4VSolid* BS_ZZ;
  G4LogicalVolume* BLV_ZZ;
  G4VSolid* CS_ZZ;
  G4LogicalVolume* CLV_ZZ;

  int k=0;

  while (barycenterListTemp)
  {
    k++;

    //3 spheres to Base, Sugar, Phosphate)
    AS_ZZ = new G4Orb("Sphere", 1.*angstrom);
    ALV_ZZ = new G4LogicalVolume(AS_ZZ,fpWaterMaterial, "ALV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Blue()));
    MyVisAtt_ZZ->SetForceSolid(true);
    ALV_ZZ->SetVisAttributes(MyVisAtt_ZZ);
    new G4PVPlacement(0,
                      G4ThreeVector(barycenterListTemp->fCenterBaseX*angstrom,
                                    barycenterListTemp->fCenterBaseY*angstrom,
                                    barycenterListTemp->fCenterBaseZ*angstrom),
                                    ALV_ZZ,
                                    "AZZ",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);

    BS_ZZ = new G4Orb("Sphere", 1.*angstrom);
    BLV_ZZ = new G4LogicalVolume(BS_ZZ,fpWaterMaterial, "BLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Red()));
    MyVisAtt_ZZ->SetForceSolid(true);
    BLV_ZZ->SetVisAttributes(MyVisAtt_ZZ);
    new G4PVPlacement(0,
                      G4ThreeVector(
                          (barycenterListTemp->fCenterPhosphateX)*angstrom,
                          (barycenterListTemp->fCenterPhosphateY)*angstrom,
                          (barycenterListTemp->fCenterPhosphateZ)*angstrom),
                          BLV_ZZ,
                          "BZZ",
                          worldLV,
                          false,
                          0,
                          fCheckOverlaps);

    CS_ZZ = new G4Orb("Sphere", 1.*angstrom);
    CLV_ZZ = new G4LogicalVolume(CS_ZZ,fpWaterMaterial, "CLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
    MyVisAtt_ZZ->SetForceSolid(true);
    CLV_ZZ->SetVisAttributes(MyVisAtt_ZZ);
    new G4PVPlacement(0,
                      G4ThreeVector(
                          barycenterListTemp->fCenterSugarX*angstrom,
                          barycenterListTemp->fCenterSugarY*angstrom,
                          barycenterListTemp->fCenterSugarZ*angstrom),
                          CLV_ZZ,
                          "CZZ",
                          worldLV,
                          false,
                          0,
                          fCheckOverlaps);

    //2 cylinders (Base<=>Sugar and Sugar<=>Phosphat)
    // Cylinder Base<=>Sugar
    tubS1_ZZ    = new G4Tubs( "Cylinder",
                              0.,
                              0.5*angstrom,
                              std::sqrt (
           (barycenterListTemp->fCenterBaseX-barycenterListTemp->fCenterSugarX)
         * (barycenterListTemp->fCenterBaseX-barycenterListTemp->fCenterSugarX)
         + (barycenterListTemp->fCenterBaseY-barycenterListTemp->fCenterSugarY)
         * (barycenterListTemp->fCenterBaseY-barycenterListTemp->fCenterSugarY)
         + (barycenterListTemp->fCenterBaseZ-barycenterListTemp->fCenterSugarZ)
         * (barycenterListTemp->fCenterBaseZ-barycenterListTemp->fCenterSugarZ)
                              ) /2 *angstrom,
                              0.,
                              2.*pi );

    tubLV1_ZZ    = new G4LogicalVolume(tubS1_ZZ,fpWaterMaterial, "tubLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Green()));
    MyVisAtt_ZZ->SetForceSolid(true);
    tubLV1_ZZ->SetVisAttributes(MyVisAtt_ZZ);

    G4double Ux= barycenterListTemp->fCenterBaseX-
        barycenterListTemp->fCenterSugarX;
    G4double Uy= barycenterListTemp->fCenterBaseY-
        barycenterListTemp->fCenterSugarY;
    G4double Uz= barycenterListTemp->fCenterBaseZ-
        barycenterListTemp->fCenterSugarZ;
    G4double llUll=std::sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

    Ux=Ux/llUll;
    Uy=Uy/llUll;
    Uz=Uz/llUll;

    G4ThreeVector direction = G4ThreeVector(Ux,Uy,Uz);
    G4double theta_euler =  direction.theta();
    G4double phi_euler   =  direction.phi();
    G4double psi_euler   = 0;

    //Warning : clhep Euler constructor build inverse matrix !
    G4RotationMatrix rotm1Inv  = G4RotationMatrix(phi_euler+pi/2,
                                                  theta_euler,
                                                  psi_euler);
    G4RotationMatrix rotm1 = rotm1Inv.inverse();
    G4ThreeVector translm1 = G4ThreeVector(
        (barycenterListTemp->fCenterBaseX+barycenterListTemp->fCenterSugarX)/2.
        *angstrom,
        (barycenterListTemp->fCenterBaseY+barycenterListTemp->fCenterSugarY)/2.
        *angstrom,
        (barycenterListTemp->fCenterBaseZ+barycenterListTemp->fCenterSugarZ)/2.
        *angstrom);
    G4Transform3D transform1 = G4Transform3D(rotm1,translm1);
    new G4PVPlacement(transform1,  // rotation translation
                      tubLV1_ZZ,"atomZZ",worldLV,false, 0, fCheckOverlaps);

    //Cylinder Sugar<=>Phosphat
    tubS2_ZZ    = new G4Tubs( "Cylinder2",
                              0.,
                              0.5*angstrom,
                              std::sqrt (
      (barycenterListTemp->fCenterSugarX-barycenterListTemp->fCenterPhosphateX)
    * (barycenterListTemp->fCenterSugarX-barycenterListTemp->fCenterPhosphateX)
    + (barycenterListTemp->fCenterSugarY-barycenterListTemp->fCenterPhosphateY)
    * (barycenterListTemp->fCenterSugarY-barycenterListTemp->fCenterPhosphateY)
    + (barycenterListTemp->fCenterSugarZ-barycenterListTemp->fCenterPhosphateZ)
    * (barycenterListTemp->fCenterSugarZ-barycenterListTemp->fCenterPhosphateZ)
                              ) /2 *angstrom,
                              0.,
                              2.*pi );

    tubLV2_ZZ    = new G4LogicalVolume(tubS2_ZZ, fpWaterMaterial, "tubLV2_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(1,0.5,0));
    MyVisAtt_ZZ->SetForceSolid(true);
    tubLV2_ZZ->SetVisAttributes(MyVisAtt_ZZ);

    Ux= barycenterListTemp->fCenterSugarX-
        barycenterListTemp->fCenterPhosphateX;
    Uy= barycenterListTemp->fCenterSugarY-
        barycenterListTemp->fCenterPhosphateY;
    Uz= barycenterListTemp->fCenterSugarZ-
        barycenterListTemp->fCenterPhosphateZ;
    llUll=std::sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

    Ux=Ux/llUll;
    Uy=Uy/llUll;
    Uz=Uz/llUll;

    direction = G4ThreeVector(Ux,Uy,Uz);
    theta_euler = direction.theta();
    phi_euler = direction.phi();
    psi_euler = 0;

    //Warning : clhep Euler constructor build inverse matrix !
    rotm1Inv  = G4RotationMatrix(phi_euler+pi/2,theta_euler,psi_euler);
    rotm1 = rotm1Inv.inverse();
    translm1 = G4ThreeVector(
        (barycenterListTemp->fCenterSugarX+
            barycenterListTemp->fCenterPhosphateX)/2.*angstrom,
            (barycenterListTemp->fCenterSugarY+
                barycenterListTemp->fCenterPhosphateY)/2.*angstrom,
                (barycenterListTemp->fCenterSugarZ+
                    barycenterListTemp->fCenterPhosphateZ)/2.*angstrom);
    transform1 = G4Transform3D(rotm1,translm1);
    new G4PVPlacement(transform1, // rotation translation
                      tubLV2_ZZ,
                      "atomZZ",
                      worldLV,
                      false,
                      0,
                      fCheckOverlaps);

    barycenterListTemp=barycenterListTemp->GetNext();
  }//end of while sur moleculeListTemp
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////
///////////////// BEGIN draw a bounding volume
//
void DetectorConstruction::DrawBoundingVolume(G4LogicalVolume* worldLV,
    Molecule *moleculeListTemp)
{
  CheckMaterials();

  double dX,dY,dZ;//Dimensions for bounding volume
  double tX,tY,tZ;//Translation for bounding volume
  fPDBlib.ComputeBoundingVolumeParams(moleculeListTemp,
                                      dX, dY, dZ,
                                      tX, tY, tZ);

  //[Geometry] Build a box :
  G4VSolid* boundingS
  = new G4Box("Bounding", dX*1*angstrom, dY*1*angstrom, dZ*1*angstrom);

  G4LogicalVolume* boundingLV
  = new G4LogicalVolume(boundingS, fpWaterMaterial, "BoundingLV");

  G4RotationMatrix *pRot = new G4RotationMatrix();

  G4VisAttributes *MyVisAtt_ZZ = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));
  boundingLV->SetVisAttributes(MyVisAtt_ZZ);

  new G4PVPlacement(pRot,
                    G4ThreeVector(
                        tX*1*angstrom,
                        tY*1*angstrom,
                        tZ*1*angstrom),  // at (x,y,z)
                        boundingLV,
                        "boundingPV",
                        worldLV
                        ,false,
                        0,
                        fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::LoadPDBfile(G4String fileName)
{
  G4cout<<"Load PDB file : "<<fileName<<"."<<G4endl<<G4endl;
  fPdbFileName=fileName;
#ifdef G4MULTITHREADED
  G4MTRunManager::GetRunManager()->DefineWorldVolume(
      DefineVolumes(fPdbFileName,fChosenOption)
  );
  G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
  G4RunManager::GetRunManager()->DefineWorldVolume(
      DefineVolumes(fPdbFileName,fChosenOption)
  );
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildBoundingVolume()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
    G4cout<<"Build only world volume and bounding volume"
        " for computation."<<G4endl<<G4endl;
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,10)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,10)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawAtoms_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,1)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,1)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawNucleotides_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,2)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,2)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawResidues_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,3)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,3)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawAtomsWithBounding_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,11)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,11)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawNucleotidesWithBounding_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,12)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,12)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawResiduesWithBounding_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,13)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,13)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}
