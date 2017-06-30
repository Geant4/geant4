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
// $Id: DicomDetectorConstruction.cc 104917 2017-06-29 07:39:14Z gcosmo $
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
#include "G4NistManager.hh"
// #include "G4SystemOfUnits.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include "DicomDetectorConstruction.hh"
#include "DicomPhantomZSliceHeader.hh"

#include "DicomRunAction.hh"
#include "DicomRun.hh"
#ifdef G4_DCMTK
#include "DicomFileMgr.hh"
#endif
#include "G4VisAttributes.hh"

using CLHEP::m;
using CLHEP::cm3;
using CLHEP::mole;
using CLHEP::g;
using CLHEP::mg;

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

   fNVoxelX(0),
   fNVoxelY(0),
   fNVoxelZ(0),
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

#ifdef G4_DCMTK
    ReadPhantomDataNew();
    ConstructPhantomContainerNew();
#else
    ReadPhantomData();
    ConstructPhantomContainer();
#endif
    
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
                                   z = 19.0, a = 39.0983* g/mole );
    G4Element* elP = new G4Element( name = "Phosphorus",
                                   symbol = "P",
                                   z = 15.0, a = 30.973976* g/mole );
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
    G4Material* trabecularBone = new G4Material( "TrabecularBone",
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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
  for( G4int ii = 0; ii < nMaterials; ii++ ){
    fin >> nmate;
    fin >> mateName;
    if( mateName[0] == '"' && mateName[mateName.length()-1] == '"' ) {
      mateName = mateName.substr(1,mateName.length()-2);
    }
    G4cout << "GmReadPhantomG4Geometry::ReadPhantomData reading nmate " << ii << " = " << nmate 
           << " mate " << mateName << G4endl;
    if( ii != nmate ) G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                                  "Wrong argument",
                                  FatalErrorInArgument,
                      "Material number should be in increasing order: wrong material number ");

    G4Material* mate = 0;
    const G4MaterialTable* matTab = G4Material::GetMaterialTable();
    std::vector<G4Material*>::const_iterator matite;
    for( matite = matTab->begin(); matite != matTab->end(); ++matite ) {
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

  fin >> fNVoxelX >> fNVoxelY >> fNVoxelZ;
  G4cout << "GmReadPhantomG4Geometry::ReadPhantomData fNVoxel X/Y/Z " << fNVoxelX << " " 
         << fNVoxelY << " " << fNVoxelZ << G4endl;
  fin >> fMinX >> fMaxX;
  fin >> fMinY >> fMaxY;
  fin >> fMinZ >> fMaxZ;
  fVoxelHalfDimX = (fMaxX-fMinX)/fNVoxelX/2.;
  fVoxelHalfDimY = (fMaxY-fMinY)/fNVoxelY/2.;
  fVoxelHalfDimZ = (fMaxZ-fMinZ)/fNVoxelZ/2.;
#ifdef G4VERBOSE
  G4cout << " Extension in X " << fMinX << " " << fMaxX << G4endl
         << " Extension in Y " << fMinY << " " << fMaxY << G4endl
         << " Extension in Z " << fMinZ << " " << fMaxZ << G4endl;
#endif

  fMateIDs = new size_t[fNVoxelX*fNVoxelY*fNVoxelZ];
  for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
    for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
      for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
        G4int mateID;
        fin >> mateID; 
        G4int nnew = ix + (iy)*fNVoxelX + (iz)*fNVoxelX*fNVoxelY;
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


void DicomDetectorConstruction::ReadVoxelDensities( std::ifstream& fin )
{
  G4String stemp;
  std::map<G4int, std::pair<G4double,G4double> > densiMinMax;
  std::map<G4int, std::pair<G4double,G4double> >::iterator mpite;
  for( size_t ii = 0; ii < thePhantomMaterialsOriginal.size(); ii++ ){
    densiMinMax[ii] = std::pair<G4double,G4double>(DBL_MAX,-DBL_MAX);
  }

  char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
  G4double densityDiff = -1.;
  if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);

  std::map<G4int,G4double> densityDiffs;
  for( size_t ii = 0; ii < thePhantomMaterialsOriginal.size(); ii++ ){
    densityDiffs[ii] = densityDiff; //currently all materials with same step
  }
  //  densityDiffs[0] = 0.0001; //air

  //--- Calculate the average material density for each material/density bin
  std::map< std::pair<G4Material*,G4int>, matInfo* > newMateDens;

  //---- Read the material densities
  G4double dens;
  for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
    for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
      for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
        fin >> dens; 
        //        G4cout << ix << " " << iy << " " << iz << " density " << dens << G4endl;
        G4int copyNo = ix + (iy)*fNVoxelX + (iz)*fNVoxelX*fNVoxelY;

        if( densityDiff != -1. ) continue; 

        //--- store the minimum and maximum density for each material (just for printing)
        mpite = densiMinMax.find( fMateIDs[copyNo] );
        if( dens < (*mpite).second.first ) (*mpite).second.first = dens;
        if( dens > (*mpite).second.second ) (*mpite).second.second = dens;
        //--- Get material from original list of material in file
        int mateID = fMateIDs[copyNo];
        std::map<G4int,G4Material*>::const_iterator imite = 
         thePhantomMaterialsOriginal.find(mateID);
        //        G4cout << copyNo << " mateID " << mateID << G4endl;
        //--- Check if density is equal to the original material density
        if( std::fabs(dens - (*imite).second->GetDensity()/CLHEP::g*CLHEP::cm3 ) < 1.e-9 ) continue;
        
        //--- Build material name with thePhantomMaterialsOriginal name + density
        //        float densityBin = densityDiffs[mateID] * (G4int(dens/densityDiffs[mateID])+0.5);
        G4int densityBin = (G4int(dens/densityDiffs[mateID]));

        G4String mateName = (*imite).second->GetName()+G4UIcommand::ConvertToString(densityBin);
        //--- Look if it is the first voxel with this material/densityBin
        std::pair<G4Material*,G4int> matdens((*imite).second, densityBin );

        std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite = 
         newMateDens.find( matdens );
        if( mppite != newMateDens.end() ){
          matInfo* mi = (*mppite).second;
          mi->fSumdens += dens;
          mi->fNvoxels++;
          fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
        } else {
          matInfo* mi = new matInfo;
          mi->fSumdens = dens;
          mi->fNvoxels = 1;
          mi->fId = newMateDens.size()+1;
          newMateDens[matdens] = mi;
          fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
        }
      }
    }
  }

  if( densityDiff != -1. ) { 
    for( mpite = densiMinMax.begin(); mpite != densiMinMax.end(); mpite++ ){
#ifdef G4VERBOSE
      G4cout << "DicomDetectorConstruction::ReadVoxelDensities ORIG MATERIALS DENSITY " 
             << (*mpite).first << " MIN " << (*mpite).second.first << " MAX " 
             << (*mpite).second.second << G4endl;
#endif
    }
  }

  //----- Build the list of phantom materials that go to Parameterisation
  //--- Add original materials
  std::map<G4int,G4Material*>::const_iterator mimite;
  for( mimite = thePhantomMaterialsOriginal.begin(); mimite != thePhantomMaterialsOriginal.end(); 
   mimite++ ){
    fMaterials.push_back( (*mimite).second );
  }
  // 
  //---- Build and add new materials
std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite;
  for( mppite= newMateDens.begin(); mppite != newMateDens.end(); mppite++ ){
    G4double averdens = (*mppite).second->fSumdens/(*mppite).second->fNvoxels;
    G4double saverdens = G4int(1000.001*averdens)/1000.;
#ifndef G4VERBOSE
    G4cout << "DicomDetectorConstruction::ReadVoxelDensities AVER DENS " << averdens << " -> " 
           << saverdens << " -> " << G4int(1000*averdens) << " " << G4int(1000*averdens)/1000 
           << " " <<  G4int(1000*averdens)/1000. << G4endl;
#endif

      G4String mateName = ((*mppite).first).first->GetName() + "_"
       + G4UIcommand::ConvertToString(saverdens);
    fMaterials.push_back( BuildMaterialWithChangingDensity( 
     (*mppite).first.first, averdens, mateName ) );
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::ReadPhantomData()
{

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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::ReadPhantomDataFile(const G4String& fname)
{
  G4cout << " DicomDetectorConstruction::ReadPhantomDataFile opening file " 
         << fname << G4endl; //GDEB
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
  char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
  G4double densityDiff = -1.;
  if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);
  if( densityDiff != -1. ) {
    for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ii++ ){
      fDensityDiffs[ii] = densityDiff; //currently all materials with 
      // same difference
    }
  }else {
    if( fMaterials.size() == 0 ) { // do it only for first slice
      for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ii++ ){
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
  G4int voxelCopyNo = (fZSliceHeaders.size()-1)*nVoxels;
  for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
    fin >> mateID;
    fMateIDs[voxelCopyNo] = mateID;
  }
  
  //----- Read material densities and build new materials if two voxels have
  //  same material but its density is in a different density interval 
  // (size of density intervals defined by densityDiff)
  G4double density;
  // number of voxels from previously read slices
  voxelCopyNo = (fZSliceHeaders.size()-1)*nVoxels;
  for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
    fin >> density;
    
    //-- Get material from list of original materials
    mateID = fMateIDs[voxelCopyNo];
    G4Material* mateOrig  = fOriginalMaterials[mateID];
    
    //-- Get density bin: middle point of the bin in which the current
    // density is included
    G4String newMateName = mateOrig->GetName();
    float densityBin = 0.;
    if( densityDiff != -1.) {
      densityBin = fDensityDiffs[mateID] * 
        (G4int(density/fDensityDiffs[mateID])+0.5);
      //-- Build the new material name
      newMateName += G4UIcommand::ConvertToString(densityBin);
    }
    
    //-- Look if a material with this name is already created
    //  (because a previous voxel was already in this density bin)
    unsigned int im;
    for( im = 0; im < fMaterials.size(); im++ ){
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
  for( unsigned int ii = 1; ii < fZSliceHeaders.size(); ii++ ) {
    *fZSliceHeaderMerged += *fZSliceHeaders[ii];
  };
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* DicomDetectorConstruction::BuildMaterialWithChangingDensity(
           const G4Material* origMate, float density, G4String newMateName )
{
  //----- Copy original material, but with new density
  G4int nelem = origMate->GetNumberOfElements();
  G4Material* mate = new G4Material( newMateName, density*g/cm3, nelem,
                                     kStateUndefined, STP_Temperature );
  
  for( G4int ii = 0; ii < nelem; ii++ ){
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
  fNVoxelX = fZSliceHeaderMerged->GetNoVoxelX();
  fNVoxelY = fZSliceHeaderMerged->GetNoVoxelY();
  fNVoxelZ = fZSliceHeaderMerged->GetNoVoxelZ();
  
  fVoxelHalfDimX = fZSliceHeaderMerged->GetVoxelHalfX();
  fVoxelHalfDimY = fZSliceHeaderMerged->GetVoxelHalfY();
  fVoxelHalfDimZ = fZSliceHeaderMerged->GetVoxelHalfZ();
#ifdef G4VERBOSE
  G4cout << " fNVoxelX " << fNVoxelX << " fVoxelHalfDimX " << fVoxelHalfDimX 
         <<G4endl;
  G4cout << " fNVoxelY " << fNVoxelY << " fVoxelHalfDimY " << fVoxelHalfDimY 
         <<G4endl;
  G4cout << " fNVoxelZ " << fNVoxelZ << " fVoxelHalfDimZ " << fVoxelHalfDimZ 
         <<G4endl;
  G4cout << " totalPixels " << fNVoxelX*fNVoxelY*fNVoxelZ <<  G4endl;
#endif
  
  //----- Define the volume that contains all the voxels
  fContainer_solid = new G4Box("phantomContainer",fNVoxelX*fVoxelHalfDimX,
                               fNVoxelY*fVoxelHalfDimY,
                               fNVoxelZ*fVoxelHalfDimZ);
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
  
  //fContainer_logic->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,0.)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::ConstructPhantomContainerNew()
{
#ifdef G4_DCMTK
  //---- Extract number of voxels and voxel dimensions
#ifdef G4VERBOSE
  G4cout << " fNVoxelX " << fNVoxelX << " fVoxelHalfDimX " << fVoxelHalfDimX 
         <<G4endl;
  G4cout << " fNVoxelY " << fNVoxelY << " fVoxelHalfDimY " << fVoxelHalfDimY 
         <<G4endl;
  G4cout << " fNVoxelZ " << fNVoxelZ << " fVoxelHalfDimZ " << fVoxelHalfDimZ 
         <<G4endl;
  G4cout << " totalPixels " << fNVoxelX*fNVoxelY*fNVoxelZ <<  G4endl;
#endif
  
  //----- Define the volume that contains all the voxels
  fContainer_solid = new G4Box("phantomContainer",fNVoxelX*fVoxelHalfDimX,
                               fNVoxelY*fVoxelHalfDimY,
                               fNVoxelZ*fVoxelHalfDimZ);
  fContainer_logic =
    new G4LogicalVolume( fContainer_solid,
   //the material is not important, it will be fully filled by the voxels
                         fMaterials[0],
                         "phantomContainer",
                         0, 0, 0 );

  G4ThreeVector posCentreVoxels((fMinX+fMaxX)/2.,(fMinY+fMaxY)/2.,(fMinZ+fMaxZ)/2.);
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
  
  //fContainer_logic->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,0.)));
#endif
}

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSDoseDeposit3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomDetectorConstruction::SetScorer(G4LogicalVolume* voxel_logic)
{

#ifdef G4VERBOSE
  G4cout << "\t SET SCORER : " << voxel_logic->GetName() << G4endl;
#endif
  
  fScorers.insert(voxel_logic);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DicomDetectorConstruction::ConstructSDandField()
{
  
#ifdef G4VERBOSE
  G4cout << "\t CONSTRUCT SD AND FIELD" << G4endl;
#endif
  
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
  G4SDManager::GetSDMpointer()->AddNewDetector( MFDet );
  //G4VPrimitiveScorer* dosedep = new G4PSDoseDeposit("DoseDeposit");
  G4VPrimitiveScorer* dosedep = 
    new G4PSDoseDeposit3D("DoseDeposit", fNVoxelX, fNVoxelY, fNVoxelZ);
  MFDet->RegisterPrimitive(dosedep);
  
  for(std::set<G4LogicalVolume*>::iterator ite = fScorers.begin(); 
      ite != fScorers.end(); ++ite) {
    SetSensitiveDetector(*ite, MFDet);
  }
  
  /*if(DicomRunAction::Instance()->GetDicomRun()) {
      DicomRunAction::Instance()->GetDicomRun()->ConstructMFD(scorer_names);
      }*/
 
}
