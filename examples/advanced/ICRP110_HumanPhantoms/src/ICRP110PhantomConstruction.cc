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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
#include "ICRP110PhantomConstruction.hh"
#include "ICRP110PhantomNestedParameterisation.hh"
#include "ICRP110PhantomMaterial_Female.hh"
#include "ICRP110PhantomMaterial_Male.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include <map>
#include <cstdlib>

ICRP110PhantomConstruction::ICRP110PhantomConstruction():
   fMotherVolume(nullptr), fPhantomContainer(nullptr),
   fNVoxelX(0), fNVoxelY(0), fNVoxelZ(0), 
   fVoxelHalfDimX(0), fVoxelHalfDimY(0), fVoxelHalfDimZ(0),
   fMinX(0),fMaxX(0), fMinY(0), fMaxY(0),
   fMinZ(0), fMaxZ(0), fNoFiles(0), fNVoxels(0),
   fMateIDs(nullptr)
{
  fMessenger = new ICRP110PhantomMessenger(this);
  // the messenger allows to set the sex of the phantom
  // interactively
  fMaterial_Female = new ICRP110PhantomMaterial_Female();
  fMaterial_Male = new ICRP110PhantomMaterial_Male();
  fSex = "female"; // Female phantom is the default option
  fSection = "head"; // Head partial phantom is the default option
}

ICRP110PhantomConstruction::~ICRP110PhantomConstruction()
{
  delete fMaterial_Female;
  delete fMaterial_Male;
  delete fMessenger;
}

G4VPhysicalVolume* ICRP110PhantomConstruction::Construct()
{
  // Define Material Air
  G4double A;  // atomic mass
  G4double Z;  // atomic number
  G4double d;  // density
  
  A = 14.01*g/mole;
  auto elN = new G4Element("Nitrogen","N",Z = 7.,A);
  A = 16.00*g/mole;
  auto elO = new G4Element("Oxygen","O",Z = 8.,A);
  
  d = 0.001 *g/cm3;
  auto matAir = new G4Material("Air",d,2);
  matAir -> AddElement(elN,0.8);
  matAir -> AddElement(elO,0.2); 

 std::vector<G4Material*> pMaterials;

 if(fSex == "female"){

  fMaterial_Female -> DefineMaterials();
//----- Store materials in a vector
    pMaterials.push_back(matAir); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("teeth")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("bone"));
    pMaterials.push_back(fMaterial_Female -> GetMaterial("humeri_upper")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("humeri_lower")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("arm_lower")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("hand")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("clavicle")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("cranium")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("femora_upper")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("femora_lower")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("leg_lower")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("foot")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("mandible")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("pelvis")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("ribs")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("scapulae")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("spine_cervical")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("spine_thoratic")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("spine_lumbar")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("sacrum")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("sternum")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("hf_upper")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("hf_lower")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("med_lowerarm")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("med_lowerleg")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("cartilage")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("skin")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("blood")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("muscle")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("liver")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("pancreas")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("brain"));
    pMaterials.push_back(fMaterial_Female -> GetMaterial("heart")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("eye")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("kidney")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("stomach")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("intestine_sml")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("intestine_lrg")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("spleen")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("thyroid")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("bladder")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("ovaries_testes")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("adrenals")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("oesophagus")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("misc")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("uterus_prostate"));
    pMaterials.push_back(fMaterial_Female -> GetMaterial("lymph")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("breast_glandular")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("breast_adipose")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("lung")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("gastro_content")); 
    pMaterials.push_back(fMaterial_Female -> GetMaterial("urine"));     
 }
 else if (fSex == "male"){
 // MATT do the same here
    fMaterial_Male -> DefineMaterials();

//----- Store materials in a vector
    pMaterials.push_back(matAir); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("teeth")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("bone"));
    pMaterials.push_back(fMaterial_Male -> GetMaterial("humeri_upper")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("humeri_lower")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("arm_lower")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("hand")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("clavicle")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("cranium")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("femora_upper")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("femora_lower")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("leg_lower")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("foot")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("mandible")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("pelvis")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("ribs")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("scapulae")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("spine_cervical")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("spine_thoratic")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("spine_lumbar")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("sacrum")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("sternum")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("hf_upper")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("hf_lower")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("med_lowerarm")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("med_lowerleg")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("cartilage")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("skin")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("blood")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("muscle")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("liver")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("pancreas")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("brain"));
    pMaterials.push_back(fMaterial_Male -> GetMaterial("heart")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("eye")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("kidney")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("stomach")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("intestine_sml")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("intestine_lrg")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("spleen")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("thyroid")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("bladder")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("ovaries_testes")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("adrenals")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("oesophagus")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("misc")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("uterus_prostate"));
    pMaterials.push_back(fMaterial_Male -> GetMaterial("lymph")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("breast_glandular")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("breast_adipose")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("lung")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("gastro_content")); 
    pMaterials.push_back(fMaterial_Male -> GetMaterial("urine")); 
 
}
  
  // World Volume
  G4double worldSize = 2.*m ;
  G4Box* world = new G4Box("world", worldSize, worldSize, worldSize);

  auto logicWorld = new G4LogicalVolume(world,
			                 matAir,
				         "logicalWorld", nullptr, nullptr,nullptr);

  fMotherVolume = new G4PVPlacement(nullptr,G4ThreeVector(),
				    "physicalWorld",
				    logicWorld,
				    nullptr,
				    false,
				    0);

  logicWorld -> SetVisAttributes(G4VisAttributes::GetInvisible());
 
  G4cout << "World has been built" << G4endl; 

  G4cout << "Phantom Sex: " << fSex << G4endl;
  G4cout << "Phantom Section: " << fSection << G4endl;
  ReadPhantomData(fSex, fSection); 
  
  G4cout << "Number of X,Y,Z voxels = " << fNVoxelX << ", " << fNVoxelY << ", " << fNVoxelZ << G4endl;
  
//----- Define the volume that contains all the voxels
  G4Box* fContainer_solid = new G4Box("phantomContainer",fNVoxelX*fVoxelHalfDimX*mm,
                               fNVoxelY*fVoxelHalfDimY*mm,
                               fNVoxelZ*fVoxelHalfDimZ*mm);
 
  auto fContainer_logic = new G4LogicalVolume( fContainer_solid,
                                                            matAir,
                                                            "phantomContainer",
                                                             nullptr, nullptr, nullptr);                                                        
    
  fMaxX = fNVoxelX*fVoxelHalfDimX*mm; // Max X along X axis of the voxelised geometry 
  fMaxY = fNVoxelY*fVoxelHalfDimY*mm; // Max Y
  fMaxZ = fNVoxelZ*fVoxelHalfDimZ*mm; // Max Z

  fMinX = -fNVoxelX*fVoxelHalfDimX*mm;// Min X 
  fMinY = -fNVoxelY*fVoxelHalfDimY*mm;// Min Y
  fMinZ = -fNVoxelZ*fVoxelHalfDimZ*mm;// Min Z

  G4ThreeVector posCentreVoxels((fMinX+fMaxX)/2.,(fMinY+fMaxY)/2.,(fMinZ+fMaxZ)/2.);

  G4cout << " placing voxel container volume at " << posCentreVoxels << G4endl;

   
  fPhantomContainer
  = new G4PVPlacement(nullptr,                     // rotation
                      posCentreVoxels,
                      fContainer_logic,     // The logic volume
                      "phantomContainer",  // Name
                      logicWorld,         // Mother
                      false,            // No op. bool.
                      1);              // Copy number
  
  fContainer_logic -> SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,0.,0.)));


// Define the voxelised phantom here
// Replication of air Phantom Volume.

//--- Slice the phantom along Y axis
   G4String yRepName("RepY");
   G4VSolid* solYRep = new G4Box(yRepName,fNVoxelX*fVoxelHalfDimX,
                                  fVoxelHalfDimY, fNVoxelZ*fVoxelHalfDimZ);
   auto logYRep = new G4LogicalVolume(solYRep,matAir,yRepName);
   new G4PVReplica(yRepName,logYRep,fContainer_logic,kYAxis, fNVoxelY,fVoxelHalfDimY*2.);
	
   logYRep -> SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));   

//--- Slice the phantom along X axis 
   G4String xRepName("RepX");
   G4VSolid* solXRep = new G4Box(xRepName,fVoxelHalfDimX,fVoxelHalfDimY,
                                  fNVoxelZ*fVoxelHalfDimZ);
   auto logXRep = new G4LogicalVolume(solXRep,matAir,xRepName);
   new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,fNVoxelX,fVoxelHalfDimX*2.);

   logXRep -> SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
	
   //----- Voxel solid and logical volumes
   //--- Slice along Z axis 
   G4VSolid* solidVoxel = new G4Box("phantom",fVoxelHalfDimX, fVoxelHalfDimY,fVoxelHalfDimZ);
   auto logicVoxel = new G4LogicalVolume(solidVoxel,matAir,"phantom");

   logicVoxel -> SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

   // Parameterisation to define the material of each voxel
   G4ThreeVector halfVoxelSize(fVoxelHalfDimX,fVoxelHalfDimY,fVoxelHalfDimZ);
      
   auto param =  new ICRP110PhantomNestedParameterisation(halfVoxelSize, pMaterials);

   new G4PVParameterised("phantom",    // their name
                          logicVoxel, // their logical volume
                          logXRep,      // Mother logical volume
                          kZAxis,       // Are placed along this axis
                          fNVoxelZ,      // Number of cells
                          param);       // Parameterisation

    param -> SetMaterialIndices(fMateIDs); // fMateIDs is  the vector with Material ID associated to each voxel, from ASCII input data files.
    param -> SetNoVoxel(fNVoxelX,fNVoxelY,fNVoxelZ);

  return fMotherVolume;
}

void ICRP110PhantomConstruction::ReadPhantomData(const G4String& sex, const G4String& section)
{

  // This method reads the information of ICRPdata/FemaleData.dat or
  // ICRPdata/MaleData.data depending on the sex of the chosen phantom

fSex = sex;
fSection = section;

G4String dataFile;

    if (fSex == "female")
    {
        if (fSection == "head")
        {
          dataFile = "ICRPdata/FemaleHead.dat";
        }
        else if (fSection == "trunk")
        {
          dataFile = "ICRPdata/FemaleTrunk.dat";
        }
        else if (fSection == "full")
        {
          dataFile = "ICRPdata/FemaleData.dat"; 
        }
    }
    if (fSex == "male")
    {
        if (fSection == "head")
        {
          dataFile = "ICRPdata/MaleHead.dat";
        }
        else if (fSection == "trunk")
        {
          dataFile = "ICRPdata/MaleTrunk.dat";
        }
        else if (fSection == "full")
        {
          dataFile = "ICRPdata/MaleData.dat"; 
        }
    }
    
    G4cout << "Data file " << dataFile << " is read by Detector Construction." << G4endl;
    
// The data.dat file in directory/build/ICRPdata/ contains the information 
// to build the phantoms. For more details look in the README file.

  //input file named finDF which consists of dataFile as a string object
  std::ifstream finDF(dataFile.c_str()); 
  
 
  G4String fname;

if(finDF.good() != 1 ) //check that the file is good and working
 { 
  G4String descript = "Problem reading data file: "+dataFile;
  G4Exception(" HumanPhantomConstruction::ReadPhantomData"," ", 
              FatalException,descript);
  }

  finDF >> fNoFiles;
  G4cout << "Number of files = " << fNoFiles << G4endl;
  finDF >> fNVoxelX;      //Inputs number of X-Voxels
  finDF >> fNVoxelY;      //Y-Voxels
  fNVoxelZ = fNoFiles;    //Z-Voxels (equal to number of slice files built/read)
  finDF >> fVoxelHalfDimX;
  finDF >> fVoxelHalfDimY;
  finDF >> fVoxelHalfDimZ;
  G4cout << "Number of X,Y,Z voxels = " << fNVoxelX << ", " << fNVoxelY << ", " << fNVoxelZ <<G4endl;
  
  fNVoxels = fNVoxelX*fNVoxelY*fNVoxelZ; 
  G4cout << "Total Number of Voxels = " << fNVoxels << G4endl;
  
  G4int nMaterials;
  finDF >> nMaterials;
  G4String mateName;
  G4int nmate;

//-----Read materials and associate with material ID number------//
  
  for( G4int ii = 0; ii < nMaterials; ii++ ){
    finDF >> nmate;
    finDF >> mateName;
    
    // This allows to skip empty spaces and tabs in the string 
    if( mateName[0] == '"' && mateName[mateName.length()-1] == '"' ) 
    {
      mateName = mateName.substr(1,mateName.length()-2); 
    }
 
    // To uncomment for eventual debugging
    /* G4cout << "GmReadPhantomG4Geometry::ReadPhantomData reading nmate " 
           << ii << " = " << nmate 
           << " mate " << mateName << G4endl;*/
 
    if( ii != nmate ) {
    G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                "Wrong argument",
                FatalErrorInArgument,
                "Material number should be in increasing order:wrong material number");
                }
      }
  
fMateIDs = new size_t[fNVoxels]; //Array with Material ID for each voxel

G4cout << "ICRP110PhantomConstruction::ReadPhantomDataFile is openining the following phantom files: " << G4endl;
  
for(G4int i = 0; i < fNoFiles; i++ )
  {
    finDF >> fname;
    ReadPhantomDataFile(fSex, fname, i); 
  }

finDF.close();
}

//----------------Opens phantom ASCII slice files to construct the phantom from-----------------//
  
void ICRP110PhantomConstruction::ReadPhantomDataFile(const G4String& sex, const G4String& fileName, G4int numberFile)
{
G4cout << fileName << G4endl;
         
fSex = sex;

G4String slice;

    if (fSex == "female")
    {
      slice = "ICRPdata/ICRP110_g4dat/AF/"+fileName;
    }
    if (fSex == "male")
    {
      slice = "ICRPdata/ICRP110_g4dat/AM/"+fileName;
    }  
  
  std::ifstream fin(slice.c_str(), std::ios_base::in);
  
  if( !fin.is_open() ) {
    G4Exception("HumanPhantomConstruction::ReadPhantomDataFile",
                "",
                FatalErrorInArgument,
                G4String("File not found " + fileName ).c_str());
  }

    for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
      for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
      if (ix == 0 && iy == 0)
        {
          G4int dudX,dudY,dudZ;      
          fin >> dudX >> dudY >> dudZ ;
          // Dummy method to skip the first three lines of the files
          // which are not used here
        }
     
        G4int nnew = ix + (iy)*fNVoxelX + numberFile*fNVoxelX*fNVoxelY;
        G4int OrgID;
        fin >> OrgID; 

        G4int mateID_out;

// The code below associates organ ID numbers (called here mateID) from ASCII slice
// files with material ID numbers (called here mateID_out) as defined in ICRP110PhantomMaterials
// Material and Organ IDs are associated as stated in AM_organs.dat and FM_organs.dat depending on
// the sex of the phantom (male and female, respctively)

	if (OrgID==128)
	{
	mateID_out=1;
	}
	
	else if (OrgID==13 || OrgID==16 || OrgID==19 || OrgID==22 || OrgID==24 || OrgID==26 || OrgID==28 || OrgID==31 || OrgID==34 || OrgID==37 || OrgID==39 || OrgID==41 || OrgID==43 || OrgID==45 || OrgID==47 || OrgID==49 || OrgID==51 || OrgID==53 || OrgID==55)
	{
	mateID_out=2;
	}

	else if (OrgID==14)
	{
	mateID_out=3;
	}

	else if (OrgID==17)
	{
	mateID_out=4;
	}

	else if (OrgID==20)
	{
	mateID_out=5;
	}

	else if (OrgID==23)
	{
	mateID_out=6;
	}

	else if (OrgID==25)
	{
	mateID_out=7;
	}

	else if (OrgID==27)
	{
	mateID_out=8;
	}

	else if (OrgID==29)
	{
	mateID_out=9;
	}

	else if (OrgID==32)
	{
	mateID_out=10;
	}

	else if (OrgID==35)
	{
	mateID_out=11;
	}
	 
	else if (OrgID==38)
	{
	mateID_out=12;
	}

	else if (OrgID==40)
	{
	mateID_out=13;
	}

	else if (OrgID==42)
	{
	mateID_out=14;
	}

	else if (OrgID==44)
	{	
	mateID_out=15;
	}

	else if (OrgID==46)
	{
	mateID_out=16;
	}

	else if (OrgID==48)
	{
	mateID_out=17;
	}

	else if (OrgID==50)
	{
	mateID_out=18;
	}

	else if (OrgID==52)
	{
	mateID_out=19;
	}

	else if (OrgID==54)
	{
	mateID_out=20;
	}

	else if (OrgID==56)
	{
	mateID_out=21;
	}

	else if (OrgID==15 || OrgID==30)
	{
	mateID_out=22;
	}

	else if (OrgID==18 || OrgID==33)
	{
	mateID_out=23;
	}

	else if (OrgID==21)
	{
	mateID_out=24;
	}

	else if (OrgID==36)
	{
	mateID_out=25;
	}

	else if (OrgID==57 || OrgID==58 || OrgID==59 || OrgID==60)	
	{
	mateID_out=26;
	}

	else if (OrgID==122 || OrgID==123 || OrgID==124 || OrgID==125 || OrgID==141 )		
	{
	mateID_out=27;
	}

	else if (OrgID==9 || OrgID==10 || OrgID==11 || OrgID==12 || OrgID==88 || OrgID==96 || OrgID==98)
	{
	mateID_out=28;
	}

	else if (OrgID==5 || OrgID==6 || OrgID==106 || OrgID==107 || OrgID==108 || OrgID==109 || OrgID==133)
	{
	mateID_out=29;
	}

	else if (OrgID==95)
	{
	mateID_out=30;
	}

	else if (OrgID==113)
	{
	mateID_out=31;
	}

	else if (OrgID==61)
	{
	mateID_out=32;
	}

	else if (OrgID==87)
	{
	mateID_out=33;
	}

	else if (OrgID==66 || OrgID==67 || OrgID==68 || OrgID==69)
	{
	mateID_out=34;
	}

	else if (OrgID==89 || OrgID==90 || OrgID==91 || OrgID==92 || OrgID==93 || OrgID==94)
	{
	mateID_out=35;
	}

	else if (OrgID==72)
	{
	mateID_out=36;
	}

	else if (OrgID==74)
	{
	mateID_out=37;
	}

	else if (OrgID==76 || OrgID==78 || OrgID==80 || OrgID==82 || OrgID==84 || OrgID==86)
	{
	mateID_out=38;
	}

	else if (OrgID==127)
	{
	mateID_out=39;
	}

	else if (OrgID==132)
	{
	mateID_out=40;
	}

	else if (OrgID==137)
	{
	mateID_out=41;
	}

	else if (OrgID==111 || OrgID==112 || OrgID==129 || OrgID==130)
	{
	mateID_out=42;
	}

	else if (OrgID==1 || OrgID==2)
	{
	mateID_out=43;
	}

	else if (OrgID==110)
	{
	mateID_out=44;
	}

	else if (OrgID==3 || OrgID==4 || OrgID==7 || OrgID==8 || OrgID==70 || OrgID==71 || OrgID==114 || OrgID==120 || OrgID==121 || OrgID==126 || OrgID==131 || OrgID==134 || OrgID==135 || OrgID == 136)
  {
  mateID_out=45;
  }

	else if (OrgID==115 || OrgID==139)
	{
	mateID_out=46;
	}

	else if (OrgID==100 || OrgID==101 || OrgID==102 || OrgID==103 || OrgID==104 || OrgID==105)
	{
	mateID_out=47;
	}

	else if (OrgID==63 || OrgID==65)
	{
	mateID_out=48;
	}

	else if (OrgID==62 || OrgID==64 || OrgID==116 || OrgID==117 || OrgID==118 || OrgID==119)
	{
	mateID_out=49;
	}

	else if (OrgID==97 || OrgID==99)
	{
	mateID_out=50;
	}

	else if (OrgID==73 || OrgID==75 || OrgID==77 || OrgID==79 || OrgID==81 || OrgID==83 || OrgID==85)
	{
	mateID_out=51;
	}

	else if (OrgID==138)
	{
	mateID_out=52;
	}

	else if (OrgID==0 || OrgID==140)
	{
	mateID_out=0;
	}

	else 
	{
	mateID_out=OrgID;
	}
        
        G4int nMaterials = 53;
        if( mateID_out < 0 || mateID_out >= nMaterials ) {
          G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                      "Wrong index in phantom file",
                      FatalException,
                      G4String("It should be between 0 and "
                              + G4UIcommand::ConvertToString(nMaterials-1) 
                              + ", while it is " 
                              + G4UIcommand::ConvertToString(OrgID)).c_str());
        
//-------------Store Material IDs and position/reference number within phantom in vector---------------//
	 }
	
          fMateIDs[nnew] = mateID_out;
	
         
      }
   }
}

//-----------Define phantom sex (male or female) to be constructed-------------//
void ICRP110PhantomConstruction::SetPhantomSex(G4String newSex)
{
  fSex = newSex;
  
  if (fSex == "male")
    {
      G4cout << ">> Male Phantom will be built." << G4endl;
    }
  if (fSex == "female")
    {
      G4cout << ">> Female Phantom will be built." << G4endl;
    }
  if ((fSex != "female") && (fSex != "male"))
    G4cout << fSex << " is not defined!" << G4endl;
} 
  
//-----------Define phantom section to be constructed-------------//
void ICRP110PhantomConstruction::SetPhantomSection(G4String newSection)
{
  fSection = newSection;
  if (fSection == "head")
    {
      G4cout << ">> Partial Head Phantom will be built." << G4endl;
    }
  if (fSection == "trunk")
    {
      G4cout << ">> Partial Trunk Phantom will be built." << G4endl;
    }
  if (fSection == "full")
    {
      G4cout << ">> Full/Custom Phantom will be built." << G4endl;
    }
  if ((fSection != "head") && (fSection != "trunk") && (fSection != "full"))
    G4cout << fSection << " is not defined!" << G4endl;  

}
