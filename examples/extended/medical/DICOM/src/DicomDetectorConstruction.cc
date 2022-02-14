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
<<<<<<< HEAD
// $Id: DicomDetectorConstruction.cc 92820 2015-09-17 15:22:14Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file  medical/DICOM/src/DicomDetectorConstruction.cc
/// \brief Implementation of the DicomDetectorConstruction class
//

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4UIcommand.hh"
#include "G4PhysicalConstants.hh"
<<<<<<< HEAD
#include "G4SystemOfUnits.hh"
=======
#include "G4NistManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#include "DicomDetectorConstruction.hh"
#include "DicomPhantomZSliceHeader.hh"
#include "DicomHandler.hh"

<<<<<<< HEAD
#include "DicomRunAction.hh"
#include "DicomRun.hh"
=======
#ifdef G4_DCMTK
#include "DicomFileMgr.hh"
#endif
#include "G4VisAttributes.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomDetectorConstruction::DicomDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAir(0),

   fWorld_solid(0),
   fWorld_logic(0),
   fWorld_phys(0),

   fContainer_solid(0),
   fContainer_logic(0),
   fContainer_phys(0),

   fNoFiles(0),
   fMateIDs(0),

   fZSliceHeaderMerged(0),

   fNoVoxelsX(0),
   fNoVoxelsY(0),
   fNoVoxelsZ(0),
   fVoxelHalfDimX(0),
   fVoxelHalfDimY(0),
   fVoxelHalfDimZ(0),

   fConstructed(false)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomDetectorConstruction::~DicomDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DicomDetectorConstruction::Construct()
{
  if(!fConstructed || fWorld_phys == 0) {
    fConstructed = true;
    InitialisationOfMaterials();

    //----- Build world
    G4double worldXDimension = 1.*m;
    G4double worldYDimension = 1.*m;
    G4double worldZDimension = 1.*m;

    fWorld_solid = new G4Box( "WorldSolid",
                             worldXDimension,
                             worldYDimension,
                             worldZDimension );

    fWorld_logic = new G4LogicalVolume( fWorld_solid,
                                       fAir,
                                       "WorldLogical",
                                       0, 0, 0 );

    fWorld_phys = new G4PVPlacement( 0,
                                    G4ThreeVector(0,0,0),
                                    "World",
                                    fWorld_logic,
                                    0,
                                    false,
                                    0 );

    ReadPhantomData();

    ConstructPhantomContainer();
    ConstructPhantom();
  }
    return fWorld_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::InitialisationOfMaterials()
{
    // Creating elements :
    G4double z, a, density;
    G4String name, symbol;

    G4Element* elC = new G4Element( name = "Carbon",
                                   symbol = "C",
                                   z = 6.0, a = 12.011 * g/mole );
    G4Element* elH = new G4Element( name = "Hydrogen",
                                   symbol = "H",
                                   z = 1.0, a = 1.008  * g/mole );
    G4Element* elN = new G4Element( name = "Nitrogen",
                                   symbol = "N",
                                   z = 7.0, a = 14.007 * g/mole );
    G4Element* elO = new G4Element( name = "Oxygen",
                                   symbol = "O",
                                   z = 8.0, a = 16.00  * g/mole );
    G4Element* elNa = new G4Element( name = "Sodium",
                                    symbol = "Na",
                                    z= 11.0, a = 22.98977* g/mole );
    G4Element* elS = new G4Element( name = "Sulfur",
                                   symbol = "S",
                                   z = 16.0,a = 32.065* g/mole );
    G4Element* elCl = new G4Element( name = "Chlorine",
                                    symbol = "P",
                                    z = 17.0, a = 35.453* g/mole );
    G4Element* elK = new G4Element( name = "Potassium",
                                   symbol = "P",
                                   z = 19.0, a = 30.0983* g/mole );
    G4Element* elP = new G4Element( name = "Phosphorus",
                                   symbol = "P",
                                   z = 30.0, a = 30.973976* g/mole );
    G4Element* elFe = new G4Element( name = "Iron",
                                    symbol = "Fe",
                                    z = 26, a = 56.845* g/mole );
    G4Element* elMg = new G4Element( name = "Magnesium",
                                    symbol = "Mg",
                                    z = 12.0, a = 24.3050* g/mole );
    G4Element* elCa = new G4Element( name="Calcium",
                                    symbol = "Ca",
                                    z = 20.0, a = 40.078* g/mole );

    // Creating Materials :
    G4int numberofElements;

    // Air
    fAir = new G4Material( "Air",
                          1.290*mg/cm3,
                          numberofElements = 2 );
    fAir->AddElement(elN, 0.7);
    fAir->AddElement(elO, 0.3);

    //  Lung Inhale
    G4Material* lunginhale = new G4Material( "LungInhale",
                                            density = 0.217*g/cm3,
                                            numberofElements = 9);
    lunginhale->AddElement(elH,0.103);
    lunginhale->AddElement(elC,0.105);
    lunginhale->AddElement(elN,0.031);
    lunginhale->AddElement(elO,0.749);
    lunginhale->AddElement(elNa,0.002);
    lunginhale->AddElement(elP,0.002);
    lunginhale->AddElement(elS,0.003);
    lunginhale->AddElement(elCl,0.002);
    lunginhale->AddElement(elK,0.003);

    // Lung exhale
    G4Material* lungexhale = new G4Material( "LungExhale",
                                            density = 0.508*g/cm3,
                                            numberofElements = 9 );
    lungexhale->AddElement(elH,0.103);
    lungexhale->AddElement(elC,0.105);
    lungexhale->AddElement(elN,0.031);
    lungexhale->AddElement(elO,0.749);
    lungexhale->AddElement(elNa,0.002);
    lungexhale->AddElement(elP,0.002);
    lungexhale->AddElement(elS,0.003);
    lungexhale->AddElement(elCl,0.002);
    lungexhale->AddElement(elK,0.003);

    // Adipose tissue
    G4Material* adiposeTissue = new G4Material( "AdiposeTissue",
                                               density = 0.967*g/cm3,
                                               numberofElements = 7);
    adiposeTissue->AddElement(elH,0.114);
    adiposeTissue->AddElement(elC,0.598);
    adiposeTissue->AddElement(elN,0.007);
    adiposeTissue->AddElement(elO,0.278);
    adiposeTissue->AddElement(elNa,0.001);
    adiposeTissue->AddElement(elS,0.001);
    adiposeTissue->AddElement(elCl,0.001);

    // Breast
    G4Material* breast = new G4Material( "Breast",
                                        density = 0.990*g/cm3,
                                        numberofElements = 8 );
    breast->AddElement(elH,0.109);
    breast->AddElement(elC,0.506);
    breast->AddElement(elN,0.023);
    breast->AddElement(elO,0.358);
    breast->AddElement(elNa,0.001);
    breast->AddElement(elP,0.001);
    breast->AddElement(elS,0.001);
    breast->AddElement(elCl,0.001);

    // Water
    G4Material* water = new G4Material( "Water",
                                       density = 1.0*g/cm3,
                                       numberofElements = 2 );
    water->AddElement(elH,0.112);
    water->AddElement(elO,0.888);

    // Muscle
    G4Material* muscle = new G4Material( "Muscle",
                                        density = 1.061*g/cm3,
                                        numberofElements = 9 );
    muscle->AddElement(elH,0.102);
    muscle->AddElement(elC,0.143);
    muscle->AddElement(elN,0.034);
    muscle->AddElement(elO,0.710);
    muscle->AddElement(elNa,0.001);
    muscle->AddElement(elP,0.002);
    muscle->AddElement(elS,0.003);
    muscle->AddElement(elCl,0.001);
    muscle->AddElement(elK,0.004);

    // Liver
    G4Material* liver = new G4Material( "Liver",
                                       density = 1.071*g/cm3,
                                       numberofElements = 9);
    liver->AddElement(elH,0.102);
    liver->AddElement(elC,0.139);
    liver->AddElement(elN,0.030);
    liver->AddElement(elO,0.716);
    liver->AddElement(elNa,0.002);
    liver->AddElement(elP,0.003);
    liver->AddElement(elS,0.003);
    liver->AddElement(elCl,0.002);
    liver->AddElement(elK,0.003);

    // Trabecular Bone
<<<<<<< HEAD
    G4Material* trabecularBone = new G4Material( "TrabecularBone",
=======
    G4Material* trabecularBone = new G4Material("TrabecularBone",
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
                                                density = 1.159*g/cm3,
                                                numberofElements = 12 );
    trabecularBone->AddElement(elH,0.085);
    trabecularBone->AddElement(elC,0.404);
    trabecularBone->AddElement(elN,0.058);
    trabecularBone->AddElement(elO,0.367);
    trabecularBone->AddElement(elNa,0.001);
    trabecularBone->AddElement(elMg,0.001);
    trabecularBone->AddElement(elP,0.034);
    trabecularBone->AddElement(elS,0.002);
    trabecularBone->AddElement(elCl,0.002);
    trabecularBone->AddElement(elK,0.001);
    trabecularBone->AddElement(elCa,0.044);
    trabecularBone->AddElement(elFe,0.001);

    // Dense Bone
    G4Material* denseBone = new G4Material( "DenseBone",
                                           density = 1.575*g/cm3,
                                           numberofElements = 11 );
    denseBone->AddElement(elH,0.056);
    denseBone->AddElement(elC,0.235);
    denseBone->AddElement(elN,0.050);
    denseBone->AddElement(elO,0.434);
    denseBone->AddElement(elNa,0.001);
    denseBone->AddElement(elMg,0.001);
    denseBone->AddElement(elP,0.072);
    denseBone->AddElement(elS,0.003);
    denseBone->AddElement(elCl,0.001);
    denseBone->AddElement(elK,0.001);
    denseBone->AddElement(elCa,0.146);

<<<<<<< HEAD
    //----- Put the materials in a vector
    fOriginalMaterials.push_back(fAir); // rho = 0.00129
    fOriginalMaterials.push_back(lunginhale); // rho = 0.217
    fOriginalMaterials.push_back(lungexhale); // rho = 0.508
    fOriginalMaterials.push_back(adiposeTissue); // rho = 0.967
    fOriginalMaterials.push_back(breast ); // rho = 0.990
    fOriginalMaterials.push_back(water); // rho = 1.018
    fOriginalMaterials.push_back(muscle); // rho = 1.061
    fOriginalMaterials.push_back(liver); // rho = 1.071
    fOriginalMaterials.push_back(trabecularBone); // rho = 1.159
    fOriginalMaterials.push_back(denseBone); // rho = 1.575

=======
    // Cortical Bone (ICRP - NIST)
    G4Material* corticalBone = new G4Material ("CorticalBone", 1.85 * g/cm3, 
                                               numberofElements = 9);
    corticalBone->AddElement(elH, 4.7234*perCent);
    corticalBone->AddElement(elC, 14.4330*perCent);
    corticalBone->AddElement(elN, 4.199*perCent);
    corticalBone->AddElement(elO, 44.6096*perCent);
    corticalBone->AddElement(elMg, 0.22*perCent);
    corticalBone->AddElement(elP, 10.497*perCent);
    corticalBone->AddElement(elS, 0.315*perCent);
    corticalBone->AddElement(elCa, 20.993*perCent);
    corticalBone->AddElement(elZn, 0.01*perCent);


    // Tooth enamel 
    G4Material* toothEnamel = new G4Material ("ToothEnamel", 2.89 * g/cm3,
                                              numberofElements = 10);
    toothEnamel->AddElement(elH, 0.95*perCent);
    toothEnamel->AddElement(elC, 1.11*perCent);
    toothEnamel->AddElement(elN, 0.23*perCent);
    toothEnamel->AddElement(elO,41.66*perCent);
    toothEnamel->AddElement(elNa, 0.79*perCent);
    toothEnamel->AddElement(elMg, 0.23*perCent);
    toothEnamel->AddElement(elP, 18.71*perCent);
    toothEnamel->AddElement(elCl, 0.34*perCent);
    toothEnamel->AddElement(elCa, 35.97*perCent);
    toothEnamel->AddElement(elZn, 0.02*perCent);

#ifdef DICOM_USE_HEAD
  //----- Put the materials in a vector HEAD PHANTOM
  fOriginalMaterials.push_back(fAir); //0.00129 g/cm3
  fOriginalMaterials.push_back(softTissue); // 1.055 g/cm3
  fOriginalMaterials.push_back(brainTissue); // 1.07 g/cm3
  fOriginalMaterials.push_back(spinalDisc); // 1.10 g/cm3
  fOriginalMaterials.push_back(trabecularBone_head); // 1.13 g/cm3
  fOriginalMaterials.push_back(toothDentin); // 1.66 g/cm3
  fOriginalMaterials.push_back(corticalBone);  // 1.75 g/cm3
  fOriginalMaterials.push_back(toothEnamel); // 2.04 g/cm3
 G4cout << "The materials of the DICOM Head have been used" << G4endl;
#else
  fOriginalMaterials.push_back(fAir); // rho = 0.00129
  fOriginalMaterials.push_back(lunginhale); // rho = 0.217
  fOriginalMaterials.push_back(lungexhale); // rho = 0.508
  fOriginalMaterials.push_back(adiposeTissue); // rho = 0.967
  fOriginalMaterials.push_back(breast ); // rho = 0.990
  fOriginalMaterials.push_back(water); // rho = 1.018
  fOriginalMaterials.push_back(muscle); // rho = 1.061
  fOriginalMaterials.push_back(liver); // rho = 1.071
  fOriginalMaterials.push_back(trabecularBone); // rho = 1.159 - HEAD PHANTOM
  fOriginalMaterials.push_back(denseBone); // rho = 1.575
  G4cout << "Default materials of the DICOM Extended examples have been used"
         << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomDetectorConstruction::ReadPhantomDataNew()
{
#ifdef G4_DCMTK
  G4String fileName = DicomFileMgr::GetInstance()->GetFileOutName();
  
  std::ifstream fin(fileName);
  std::vector<G4String> wl;
  G4int nMaterials;
  fin >> nMaterials;
  G4String mateName;
  G4int nmate;
  for( G4int ii = 0; ii < nMaterials; ++ii )
  {
    fin >> nmate;
    fin >> mateName;
    if( mateName[0] == '"' && mateName[mateName.length()-1] == '"' ) {
      mateName = mateName.substr(1,mateName.length()-2);
    }
    G4cout << "GmReadPhantomG4Geometry::ReadPhantomData reading nmate " 
           << ii << " = " << nmate 
           << " mate " << mateName << G4endl;
    if( ii != nmate ) 
    G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                "Wrong argument", FatalErrorInArgument,
        "Material number should be in increasing order:wrong material number");

    G4Material* mate = 0;
    const G4MaterialTable* matTab = G4Material::GetMaterialTable();
    for( auto matite = matTab->cbegin(); matite != matTab->cend(); ++matite )
    {
      if( (*matite)->GetName() == mateName ) {
        mate = *matite;
      }
    }
    if( mate == 0 ) {
      mate = G4NistManager::Instance()->FindOrBuildMaterial(mateName);
    }
    if( !mate ) G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                            "Wrong argument",
                            FatalErrorInArgument,
                            ("Material not found" + mateName).c_str());
    thePhantomMaterialsOriginal[nmate] = mate;
  }

  fin >> fNoVoxelsX >> fNoVoxelsY >> fNoVoxelsZ;
  G4cout << "GmReadPhantomG4Geometry::ReadPhantomData fNoVoxels X/Y/Z " 
         << fNoVoxelsX << " " 
         << fNoVoxelsY << " " << fNoVoxelsZ << G4endl;
  fin >> fMinX >> fMaxX;
  fin >> fMinY >> fMaxY;
  fin >> fMinZ >> fMaxZ;
  fVoxelHalfDimX = (fMaxX-fMinX)/fNoVoxelsX/2.;
  fVoxelHalfDimY = (fMaxY-fMinY)/fNoVoxelsY/2.;
  fVoxelHalfDimZ = (fMaxZ-fMinZ)/fNoVoxelsZ/2.;
#ifdef G4VERBOSE
  G4cout << " Extension in X " << fMinX << " " << fMaxX << G4endl
         << " Extension in Y " << fMinY << " " << fMaxY << G4endl
         << " Extension in Z " << fMinZ << " " << fMaxZ << G4endl;
#endif

  fMateIDs = new size_t[fNoVoxelsX*fNoVoxelsY*fNoVoxelsZ];
  for( G4int iz = 0; iz < fNoVoxelsZ; ++iz ) {
    for( G4int iy = 0; iy < fNoVoxelsY; ++iy ) {
      for( G4int ix = 0; ix < fNoVoxelsX; ++ix ) {
        G4int mateID;
        fin >> mateID; 
        G4int nnew = ix + (iy)*fNoVoxelsX + (iz)*fNoVoxelsX*fNoVoxelsY;
        if( mateID < 0 || mateID >= nMaterials ) {
          G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                      "Wrong index in phantom file",
                      FatalException,
                      G4String("It should be between 0 and "
                              + G4UIcommand::ConvertToString(nMaterials-1) 
                              + ", while it is " 
                              + G4UIcommand::ConvertToString(mateID)).c_str());
        }
        fMateIDs[nnew] = mateID;
      }
    }
  }

  ReadVoxelDensities( fin );

  fin.close();
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void DicomDetectorConstruction::ReadVoxelDensities( std::ifstream& fin )
{
  G4String stemp;
  std::map<G4int, std::pair<G4double,G4double> > densiMinMax;
  std::map<G4int, std::pair<G4double,G4double> >::iterator mpite;
  for( G4int ii = 0; ii < G4int(thePhantomMaterialsOriginal.size()); ++ii )
  {
    densiMinMax[ii] = std::pair<G4double,G4double>(DBL_MAX,-DBL_MAX);
  }

  char* part = std::getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
  G4double densityDiff = -1.;
  if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);

  std::map<G4int,G4double> densityDiffs;
  for( G4int ii = 0; ii < G4int(thePhantomMaterialsOriginal.size()); ++ii )
  {
    densityDiffs[ii] = densityDiff; //currently all materials with same step
  }
  //  densityDiffs[0] = 0.0001; //air

  //--- Calculate the average material density for each material/density bin
  std::map< std::pair<G4Material*,G4int>, matInfo* > newMateDens;

  //---- Read the material densities
  G4double dens;
  for( G4int iz = 0; iz < fNoVoxelsZ; ++iz ) {
    for( G4int iy = 0; iy < fNoVoxelsY; ++iy ) {
      for( G4int ix = 0; ix < fNoVoxelsX; ++ix ) {
        fin >> dens; 
        G4int copyNo = ix + (iy)*fNoVoxelsX + (iz)*fNoVoxelsX*fNoVoxelsY;

        if( densityDiff != -1. ) continue; 

        //--- store the minimum and maximum density for each material
        mpite = densiMinMax.find( G4int(fMateIDs[copyNo]) );
        if( dens < (*mpite).second.first ) (*mpite).second.first = dens;
        if( dens > (*mpite).second.second ) (*mpite).second.second = dens;
        //--- Get material from original list of material in file
        G4int mateID = G4int(fMateIDs[copyNo]);
        std::map<G4int,G4Material*>::const_iterator imite = 
         thePhantomMaterialsOriginal.find(mateID);

        //--- Check if density is equal to the original material density
        if(std::fabs(dens - (*imite).second->GetDensity()/CLHEP::g*CLHEP::cm3 )
           < 1.e-9 ) continue;
        
        //--- Build material name with thePhantomMaterialsOriginal name+density
        G4int densityBin = (G4int(dens/densityDiffs[mateID]));

        G4String mateName = (*imite).second->GetName()
                          + G4UIcommand::ConvertToString(densityBin);
        //--- Look if it is the first voxel with this material/densityBin
        std::pair<G4Material*,G4int> matdens((*imite).second, densityBin );

        auto mppite = newMateDens.find( matdens );
        if( mppite != newMateDens.cend() ){
          matInfo* mi = (*mppite).second;
          mi->fSumdens += dens;
          mi->fNvoxels++;
          fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
        } else {
          matInfo* mi = new matInfo;
          mi->fSumdens = dens;
          mi->fNvoxels = 1;
          mi->fId = G4int(newMateDens.size()+1);
          newMateDens[matdens] = mi;
          fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
        }
      }
    }
  }

  if( densityDiff != -1. ) { 
    for( mpite = densiMinMax.begin(); mpite != densiMinMax.end(); ++mpite )
    {
#ifdef G4VERBOSE
      G4cout << "DicomDetectorConstruction::ReadVoxelDensities"
             << " ORIG MATERIALS DENSITY " 
             << (*mpite).first << " MIN " << (*mpite).second.first << " MAX " 
             << (*mpite).second.second << G4endl;
#endif
    }
  }

  //----- Build the list of phantom materials that go to Parameterisation
  //--- Add original materials
  for( auto mimite = thePhantomMaterialsOriginal.cbegin();
       mimite != thePhantomMaterialsOriginal.cend(); ++mimite )
  {
    fMaterials.push_back( (*mimite).second );
  }
  // 
  //---- Build and add new materials
  for( auto mppite= newMateDens.cbegin(); mppite!=newMateDens.cend(); ++mppite )
  {
    G4double averdens = (*mppite).second->fSumdens/(*mppite).second->fNvoxels;
    G4double saverdens = G4int(1000.001*averdens)/1000.;
#ifndef G4VERBOSE
    G4cout << "DicomDetectorConstruction::ReadVoxelDensities AVER DENS " 
           << averdens << " -> " 
           << saverdens << " -> " << G4int(1000*averdens) << " " 
           << G4int(1000*averdens)/1000 
           << " " <<  G4int(1000*averdens)/1000. << G4endl;
#endif

    G4String mateName = ((*mppite).first).first->GetName() + "_"
       + G4UIcommand::ConvertToString(saverdens);
    fMaterials.push_back( BuildMaterialWithChangingDensity( 
     (*mppite).first.first, G4float(averdens), mateName ) );
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::ReadPhantomData()
{
<<<<<<< HEAD
=======
  G4String dataFile = DicomHandler::GetDicomDataFile();
  std::ifstream finDF(dataFile.c_str());
  G4String fname;

if(finDF.good() != 1 ) 
 { 
  G4String descript = "Problem reading data file: "+dataFile;
  G4Exception(" DicomDetectorConstruction::ReadPhantomData"," ", 
              FatalException,descript);
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

    G4String dataFile = "Data.dat";
    std::ifstream finDF(dataFile.c_str());
    G4String fname;
    if(finDF.good() != 1 ) {
        G4String descript = "Problem reading data file: "+dataFile;
        G4Exception(" DicomDetectorConstruction::ReadPhantomData",
                    "",
                    FatalException,
                    descript);
    }

<<<<<<< HEAD
    G4int compression;
    finDF >> compression; // not used here

    finDF >> fNoFiles;
    for(G4int i = 0; i < fNoFiles; i++ ) {
        finDF >> fname;
        //--- Read one data file
        fname += ".g4dcm";
        ReadPhantomDataFile(fname);
    }

    //----- Merge data headers
    MergeZSliceHeaders();

    finDF.close();
=======
  for(G4int i = 0; i < fNoFiles; ++i )
  {

    finDF >> fname;

    //--- Read one data file
    fname += ".g4dcm";

    ReadPhantomDataFile(fname);
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

  //----- Merge data headers
  MergeZSliceHeaders();
  finDF.close();
}

<<<<<<< HEAD
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
void DicomDetectorConstruction::ReadPhantomDataFile(const G4String& fname)
{
#ifdef G4VERBOSE
  G4cout << " DicomDetectorConstruction::ReadPhantomDataFile opening file " 
         << fname << G4endl;
#endif
  std::ifstream fin(fname.c_str(), std::ios_base::in);
  if( !fin.is_open() ) {
    G4Exception("DicomDetectorConstruction::ReadPhantomDataFile",
                "",
                FatalErrorInArgument,
                G4String("File not found " + fname ).c_str());
  }
  //----- Define density differences (maximum density difference to create
  // a new material)
  char* part = std::getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
  G4double densityDiff = -1.;
  if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);
  if( densityDiff != -1. )
  {
    for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ++ii )
    {
      fDensityDiffs[ii] = densityDiff; //currently all materials with 
      // same difference
    }
  }
  else
  {
    if( fMaterials.size() == 0 ) { // do it only for first slice
      for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ++ii )
      {
        fMaterials.push_back( fOriginalMaterials[ii] );
      }
    }
  }
  
  //----- Read data header
  DicomPhantomZSliceHeader* sliceHeader = new DicomPhantomZSliceHeader( fin );
  fZSliceHeaders.push_back( sliceHeader );
  
  //----- Read material indices
  G4int nVoxels = sliceHeader->GetNoVoxels();
  
  //--- If first slice, initiliaze fMateIDs
  if( fZSliceHeaders.size() == 1 ) {
    //fMateIDs = new unsigned int[fNoFiles*nVoxels];
    fMateIDs = new size_t[fNoFiles*nVoxels];
    
  }
  
  unsigned int mateID;
  // number of voxels from previously read slices
  G4int voxelCopyNo = G4int((fZSliceHeaders.size()-1)*nVoxels);
  for( G4int ii = 0; ii < nVoxels; ++ii, voxelCopyNo++ )
  {
    fin >> mateID;
    fMateIDs[voxelCopyNo] = mateID;
  }
  
  //----- Read material densities and build new materials if two voxels have
  //  same material but its density is in a different density interval 
  // (size of density intervals defined by densityDiff)
  G4double density;
  // number of voxels from previously read slices
  voxelCopyNo = G4int((fZSliceHeaders.size()-1)*nVoxels);
  for( G4int ii = 0; ii < nVoxels; ++ii, voxelCopyNo++ )
  {
    fin >> density;
    
    //-- Get material from list of original materials
    mateID = unsigned(fMateIDs[voxelCopyNo]);
    G4Material* mateOrig  = fOriginalMaterials[mateID];
    
    //-- Get density bin: middle point of the bin in which the current
    // density is included
    G4String newMateName = mateOrig->GetName();
    G4float densityBin = 0.;
    if( densityDiff != -1.) {
      densityBin = G4float(fDensityDiffs[mateID]) * 
                   (G4int(density/fDensityDiffs[mateID])+0.5);
      //-- Build the new material name
      newMateName += G4UIcommand::ConvertToString(densityBin);
    }
    
    //-- Look if a material with this name is already created
    //  (because a previous voxel was already in this density bin)
    unsigned int im;
    for( im = 0; im < fMaterials.size(); ++im )
    {
      if( fMaterials[im]->GetName() == newMateName ) {
        break;
      }
    }
    //-- If material is already created use index of this material
    if( im != fMaterials.size() ) {
      fMateIDs[voxelCopyNo] = im;
      //-- else, create the material
    } else {
      if( densityDiff != -1.) {
        fMaterials.push_back( BuildMaterialWithChangingDensity( mateOrig,
                                                  densityBin, newMateName ) );
        fMateIDs[voxelCopyNo] = fMaterials.size()-1;
      } else {
        G4cerr << " im " << im << " < " << fMaterials.size() << " name "
               << newMateName << G4endl;
        G4Exception("DicomDetectorConstruction::ReadPhantomDataFile",
                    "",
                    FatalErrorInArgument,
                    "Wrong index in material"); //it should never reach here
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::MergeZSliceHeaders()
{
  //----- Images must have the same dimension ...
  fZSliceHeaderMerged = new DicomPhantomZSliceHeader( *fZSliceHeaders[0] );
  for( unsigned int ii = 1; ii < fZSliceHeaders.size(); ++ii )
  {
    *fZSliceHeaderMerged += *fZSliceHeaders[ii];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* DicomDetectorConstruction::BuildMaterialWithChangingDensity(
           const G4Material* origMate, G4float density, G4String newMateName )
{
  //----- Copy original material, but with new density
  G4int nelem = G4int(origMate->GetNumberOfElements());
  G4Material* mate = new G4Material( newMateName, density*g/cm3, nelem,
                                     kStateUndefined, STP_Temperature );
  
  for( G4int ii = 0; ii < nelem; ++ii )
  {
    G4double frac = origMate->GetFractionVector()[ii];
    G4Element* elem = const_cast<G4Element*>(origMate->GetElement(ii));
    mate->AddElement( elem, frac );
  }
  
  return mate;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::ConstructPhantomContainer()
{
  //---- Extract number of voxels and voxel dimensions
  fNoVoxelsX = fZSliceHeaderMerged->GetNoVoxelsX();
  fNoVoxelsY = fZSliceHeaderMerged->GetNoVoxelsY();
  fNoVoxelsZ = fZSliceHeaderMerged->GetNoVoxelsZ();
  
  fVoxelHalfDimX = fZSliceHeaderMerged->GetVoxelHalfX();
  fVoxelHalfDimY = fZSliceHeaderMerged->GetVoxelHalfY();
  fVoxelHalfDimZ = fZSliceHeaderMerged->GetVoxelHalfZ();
#ifdef G4VERBOSE
  G4cout << " fNoVoxelsX " << fNoVoxelsX << " fVoxelHalfDimX " << fVoxelHalfDimX 
         <<G4endl;
  G4cout << " fNoVoxelsY " << fNoVoxelsY << " fVoxelHalfDimY " << fVoxelHalfDimY 
         <<G4endl;
  G4cout << " fNoVoxelsZ " << fNoVoxelsZ << " fVoxelHalfDimZ " << fVoxelHalfDimZ 
         <<G4endl;
  G4cout << " totalPixels " << fNoVoxelsX*fNoVoxelsY*fNoVoxelsZ <<  G4endl;
#endif
  
  //----- Define the volume that contains all the voxels
  fContainer_solid = new G4Box("phantomContainer",fNoVoxelsX*fVoxelHalfDimX,
                               fNoVoxelsY*fVoxelHalfDimY,
                               fNoVoxelsZ*fVoxelHalfDimZ);
  fContainer_logic =
    new G4LogicalVolume( fContainer_solid,
   //the material is not important, it will be fully filled by the voxels
                         fMaterials[0],
                         "phantomContainer",
                         0, 0, 0 );
  //--- Place it on the world
  G4double fOffsetX = (fZSliceHeaderMerged->GetMaxX() + 
                       fZSliceHeaderMerged->GetMinX() ) /2.;
  G4double fOffsetY = (fZSliceHeaderMerged->GetMaxY() + 
                       fZSliceHeaderMerged->GetMinY() ) /2.;
  G4double fOffsetZ = (fZSliceHeaderMerged->GetMaxZ() + 
                       fZSliceHeaderMerged->GetMinZ() ) /2.;
  G4ThreeVector posCentreVoxels(fOffsetX,fOffsetY,fOffsetZ);
#ifdef G4VERBOSE
  G4cout << " placing voxel container volume at " << posCentreVoxels << G4endl;
#endif
  fContainer_phys =
    new G4PVPlacement(0,  // rotation
                      posCentreVoxels,
                      fContainer_logic,     // The logic volume
                      "phantomContainer",  // Name
                      fWorld_logic,  // Mother
                      false,           // No op. bool.
                      1);              // Copy number
}

<<<<<<< HEAD
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::ConstructPhantomContainerNew()
{
#ifdef G4_DCMTK
  //---- Extract number of voxels and voxel dimensions
#ifdef G4VERBOSE
  G4cout << " fNoVoxelsX " << fNoVoxelsX << " fVoxelHalfDimX " << fVoxelHalfDimX 
         <<G4endl;
  G4cout << " fNoVoxelsY " << fNoVoxelsY << " fVoxelHalfDimY " << fVoxelHalfDimY 
         <<G4endl;
  G4cout << " fNoVoxelsZ " << fNoVoxelsZ << " fVoxelHalfDimZ " << fVoxelHalfDimZ 
         <<G4endl;
  G4cout << " totalPixels " << fNoVoxelsX*fNoVoxelsY*fNoVoxelsZ <<  G4endl;
#endif
  
  //----- Define the volume that contains all the voxels
  fContainer_solid = new G4Box("phantomContainer",fNoVoxelsX*fVoxelHalfDimX,
                               fNoVoxelsY*fVoxelHalfDimY,
                               fNoVoxelsZ*fVoxelHalfDimZ);
  fContainer_logic =
    new G4LogicalVolume( fContainer_solid,
   //the material is not important, it will be fully filled by the voxels
                         fMaterials[0],
                         "phantomContainer",
                         0, 0, 0 );

  G4ThreeVector posCentreVoxels((fMinX+fMaxX)/2.,
                                (fMinY+fMaxY)/2.,
                                (fMinZ+fMaxZ)/2.);
#ifdef G4VERBOSE
  G4cout << " placing voxel container volume at " << posCentreVoxels << G4endl;
#endif
  fContainer_phys =
    new G4PVPlacement(0,  // rotation
                      posCentreVoxels,
                      fContainer_logic,     // The logic volume
                      "phantomContainer",  // Name
                      fWorld_logic,  // Mother
                      false,           // No op. bool.
                      1);              // Copy number
#endif
}

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSDoseDeposit3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::SetScorer(G4LogicalVolume* voxel_logic)
{
  
  G4cout << "\n\n\n\n\t SET SCORER : " << voxel_logic->GetName() 
         << " \n\n\n" << G4endl;
  
  fScorers.insert(voxel_logic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DicomDetectorConstruction::ConstructSDandField()
{
  
  G4cout << "\n\n\n\n\t CONSTRUCT SD AND FIELD \n\n\n" << G4endl;
  
  //G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  //SDman->SetVerboseLevel(1);
  
  //
  // Sensitive Detector Name
  G4String concreteSDname = "phantomSD";
  std::vector<G4String> scorer_names;
  scorer_names.push_back(concreteSDname);
  //------------------------
  // MultiFunctionalDetector
  //------------------------
  //
  // Define MultiFunctionalDetector with name.
  // declare MFDet as a MultiFunctionalDetector scorer
  G4MultiFunctionalDetector* MFDet = 
    new G4MultiFunctionalDetector(concreteSDname);
  //SDman->AddNewDetector( MFDet );                 // Register SD to SDManager
  //G4VPrimitiveScorer* dosedep = new G4PSDoseDeposit("DoseDeposit");
  G4VPrimitiveScorer* dosedep = 
    new G4PSDoseDeposit3D("DoseDeposit", fNoVoxelsX, fNoVoxelsY, fNoVoxelsZ);
  MFDet->RegisterPrimitive(dosedep);
  
  for(auto ite = fScorers.cbegin(); ite != fScorers.cend(); ++ite)
  {
    SetSensitiveDetector(*ite, MFDet);
  }
}
