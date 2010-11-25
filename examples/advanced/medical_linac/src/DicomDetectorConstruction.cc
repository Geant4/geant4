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
// The code is part of the DICOM extended example and it was modified by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "DicomDetectorConstruction.hh"
#include "DicomPatientZSliceHeader.hh"

//-------------------------------------------------------------
DicomDetectorConstruction::DicomDetectorConstruction()
{
	this->dicomMessenger=new CML2_DicomMessenger(this);
	// added ML2 start
	this->matNames=new G4String[1000]; 
	this->nMatNames=0;
	// end
}

//-------------------------------------------------------------
DicomDetectorConstruction::~DicomDetectorConstruction()
{
}
//-------------------------------------------------------------
void DicomDetectorConstruction::Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG)
{

	this->world_phys=PVWorld;
	this->saving_in_ROG_Voxels_every_events=saving_in_ROG_Voxels_every_events;
	this->seed=seed;
	this->ROGOutFile=ROGOutFile;
	this->bSaveROG=bSaveROG;

  InitialisationOfMaterials();

  //----- Build world

  this->world_logic=this->world_phys->GetLogicalVolume();

//  G4VSolid *so= this->world_phys->GetLogicalVolume()->GetSolid();

//  G4Box* bx= (G4Box *)so;

//   G4double worldXDimension = bx->GetXHalfLength();
//   G4double worldYDimension = bx->GetYHalfLength();
//   G4double worldZDimension = bx->GetZHalfLength();

  ReadPatientData(this->dicomDirectory, this->fileData);

  ConstructPatientContainer();

    // phantom size and position
	G4ThreeVector halfSize, centre; // centre and half size of the phantom container
	halfSize.set(this->container_solid->GetXHalfLength(), this->container_solid->GetYHalfLength(), this->container_solid->GetZHalfLength());

/*    G4double offsetX = (this->fZSliceHeaderMerged->GetMaxX() + this->fZSliceHeaderMerged->GetMinX() ) /2.;
    G4double offsetY = (this->fZSliceHeaderMerged->GetMaxY() + this->fZSliceHeaderMerged->GetMinY() ) /2.;
    G4double offsetZ = (this->fZSliceHeaderMerged->GetMaxZ() + this->fZSliceHeaderMerged->GetMinZ() ) /2.;
    centre.set(offsetX,offsetY,offsetZ);*/

    centre.set(0.,0.,0.); // ML2 

std::cout <<"DICOM ROG centre " <<centre << G4endl;
std::cout <<"DICOM ROG halfSize " <<halfSize << G4endl;
std::cout <<"DICOM ROG this->nROGVoxelsX " <<this->nROGVoxelsX << G4endl;
std::cout <<"DICOM ROG this->nROGVoxelsY " <<this->nROGVoxelsY << G4endl;
std::cout <<"DICOM ROG this->nROGVoxelsZ " <<this->nROGVoxelsZ << G4endl;
	// Sensitive detector
	this->sensDet=new CML2SDWithVoxels("DicomPhantomROG", this->saving_in_ROG_Voxels_every_events, this->seed, this->ROGOutFile, this->bSaveROG, centre, halfSize, this->nROGVoxelsX, this->nROGVoxelsY, this->nROGVoxelsZ);
	G4SDManager *SDManager=G4SDManager::GetSDMpointer();
	SDManager->AddNewDetector(this->sensDet);
	
	// Read Out Geometry
	CML2ReadOutGeometry *ROG = new CML2ReadOutGeometry();
	ROG->setBuildData(this->world_phys->GetFrameTranslation(), halfSize, this->nROGVoxelsX, this->nROGVoxelsY, this->nROGVoxelsZ);
	ROG->BuildROGeometry();
	this->sensDet->SetROgeometry(ROG);

		// Region for cuts
	this->regVol= new G4Region("dicomRegion");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(.001*cm);
	this->regVol->SetProductionCuts(cuts);

  ConstructPatient();

}
// ML2 Added code start
void DicomDetectorConstruction::writeInfo()
{
	std::cout <<"\n\n\tnumber of Read Out Geometry voxels: X) "<< this->nROGVoxelsX<< "\tY) "<< this->nROGVoxelsY<<"\tZ) "<<this->nROGVoxelsZ<< G4endl;
	std::cout <<"\tnumber of dicom voxels: X) "<< this->nVoxelX<< "\tY) "<< this->nVoxelY<<"\tZ) "<<this->nVoxelZ<< G4endl;
	std::cout << "\tProduction cuts: gamma)"<<this->regVol->GetProductionCuts()->GetProductionCut("gamma")/mm<<" [mm]\t";
	std::cout << "e-)"<<this->regVol->GetProductionCuts()->GetProductionCut("e-")/mm<<" [mm]\t";
	std::cout << "e+)"<<this->regVol->GetProductionCuts()->GetProductionCut("e+")/mm<<" [mm]\n"<<G4endl;
}

G4ThreeVector DicomDetectorConstruction::getHalfContainerSize()
{
	G4double minX, minY, minZ, maxX, maxY, maxZ;
	minX=minY=minZ=maxX=maxY=maxZ=0.;
	G4int appoI1, appoI2;
	G4String appoS;
	G4double appoMin, appoMax;
	static bool bFirstTime=true;

	G4String fullDataFileName=this->dicomDirectory+this->fileData;
	std::ifstream finDF(fullDataFileName);
	G4String fname;
	if(finDF.good() != 1 ) {
	G4Exception(" DicomDetectorConstruction::ReadPatientData.  Problem reading data file: "+this->fileData);
	}
	
	G4int compression;
	finDF >> compression; // not used here
	finDF >> this->fNoFiles;
	
	G4String fullDicomFileName; 
	for(G4int i = 0; i < this->fNoFiles; i++ ) 
	{
		finDF >> fname;
		//--- Read one data file
		fname += ".g4dcm";
		fullDicomFileName=this->dicomDirectory+fname; // store the .g4dcm files name ML2
		std::ifstream fin(fullDicomFileName);
		if(fin.good() != 1 )
		{
			G4Exception(" DicomDetectorConstruction::ReadPatientData.  Problem reading data file: "+fullDicomFileName);
		}
		else
		{
			fin>>appoI1;
			for (int i=0;i<appoI1;i++)
			{
				fin >> appoI2 >> appoS;
			}
			fin >> appoI1 >> appoI1 >> appoI1;
			if (bFirstTime)
			{
				bFirstTime=false;
				fin >> appoMin >> appoMax;
				minX=appoMin;
				maxX=appoMax;
				fin >> appoMin >> appoMax;
				minY=appoMin;
				maxY=appoMax;
				fin >> appoMin >> appoMax;
				minZ=appoMin;
				maxZ=appoMax;
			}
			else
			{
				fin >> appoMin >> appoMax;
				if (minX>appoMin){minX=appoMin;}
				if (maxX<appoMax){maxX=appoMax;}
				fin >> appoMin >> appoMax;
				if (minY>appoMin){minY=appoMin;}
				if (maxY<appoMax){maxY=appoMax;}
				fin >> appoMin >> appoMax;
				if (minZ>appoMin){minZ=appoMin;}
				if (maxZ<appoMax){maxZ=appoMax;}
			}
		}
		fin.close();
	}
	this->halfContainerSize.set((maxX-minX)/2., (maxY-minY)/2., (maxZ-minZ)/2.);
	std::cout <<'\n'<<"Dicom halfContainerSize [mm]: " <<this->halfContainerSize/mm <<'\n'<< G4endl;
	return this->halfContainerSize;
}
// ML2 Added code end

//-------------------------------------------------------------
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
  air = new G4Material( "Air",
                        1.290*mg/cm3,
                        numberofElements = 2 );
  air->AddElement(elN, 0.7);
  air->AddElement(elO, 0.3); 

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
  fOriginalMaterials.push_back(air); // rho = 0.00129
  fOriginalMaterials.push_back(lunginhale); // rho = 0.217
  fOriginalMaterials.push_back(lungexhale); // rho = 0.508
  fOriginalMaterials.push_back(adiposeTissue); // rho = 0.967
  fOriginalMaterials.push_back(breast ); // rho = 0.990
  fOriginalMaterials.push_back(water); // rho = 1.018
  fOriginalMaterials.push_back(muscle); // rho = 1.061
  fOriginalMaterials.push_back(liver); // rho = 1.071
  fOriginalMaterials.push_back(trabecularBone); // rho = 1.159
  fOriginalMaterials.push_back(denseBone); // rho = 1.575

  this->ML2_G4NIST_Materials();
}

// added ML2 start

// added G4-NIST and other materials 
#include "G4NistManager.hh"

void DicomDetectorConstruction::ML2_G4NIST_Materials()
{
	G4double A, Z;
	G4int natoms;
 
	A = 1.01*g/mole;
	G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);
	A = 12.011*g/mole;
	G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);  
	A = 14.01*g/mole;
	G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);
	A = 16.00*g/mole;
	G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

	G4double dPolyurethane=1.05*g/cm3;
	G4Material* polyurethanePieno = new G4Material("polyurethanePieno", dPolyurethane,4);
	polyurethanePieno->AddElement(elH, natoms=42);
	polyurethanePieno->AddElement(elC, natoms=25);
	polyurethanePieno->AddElement(elN, natoms=2);
	polyurethanePieno->AddElement(elO, natoms=6);
    fOriginalMaterials.push_back(polyurethanePieno);

	G4Material *newNISTMaterial;
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_A-150_TISSUE"); // 1.127
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP"); // 0.92
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_B-100_BONE"); // 1.45
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_BLOOD_ICRP"); // 1.06
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU"); // 1.85
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_BONE_CORTICAL_ICRP"); // 1.85
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_BRAIN_ICRP"); // 1.03
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_EYE_LENS_ICRP"); // 1.1
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_LUNG_ICRP"); // 1.05
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_MS20_TISSUE"); // 1.0
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_MUSCLE_SKELETAL_ICRP"); // 1.04 
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_MUSCLE_STRIATED_ICRU"); // 1.04
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_MUSCLE_WITH_SUCROSE"); // 1.11
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_MUSCLE_WITHOUT_SUCROSE"); // 1.07
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_SKIN_ICRP"); // 1.1
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_SUCROSE"); // 1.5805
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP"); // 1.0
    fOriginalMaterials.push_back(newNISTMaterial);
	newNISTMaterial=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER_VAPOR"); // 0.000756182
    fOriginalMaterials.push_back(newNISTMaterial);
}

// ML2 end
//-------------------------------------------------------------
void DicomDetectorConstruction::ReadPatientData(G4String directory, G4String fileName)
{
	// changed to read and store the materials name and densities  ML2

	G4String fullDataFileName=directory+fileName;

  std::ifstream finDF(fullDataFileName);
  G4String fname;
  if(finDF.good() != 1 ) {
    G4Exception(" DicomDetectorConstruction::ReadPatientData.  Problem reading data file: "+fileName);
  }

  
  G4int compression;
  finDF >> compression; // not used here

  finDF >> this->fNoFiles;



  G4String *fullDicomFileName=new G4String[fNoFiles]; 
  for(G4int i = 0; i < fNoFiles; i++ ) {    
    finDF >> fname;
    //--- Read one data file
    fname += ".g4dcm";
	fullDicomFileName[i]=directory+fname; // store the .g4dcm files name ML2
  }

  finDF >> this->nInputMatNames;
  this->inputMatNames=new G4String[this->nInputMatNames];
  this->density=new G4double[this->nInputMatNames];

  for(G4int i = 0; i < this->nInputMatNames; i++ ) {    
	  finDF >> this->inputMatNames[i] >> this->density[i]; // store the material name and its density ML2
  }

  finDF.close();

  for(G4int i = 0; i <  this->fNoFiles; i++ ) {    
    ReadPatientDataFile(fullDicomFileName[i]);
  }
  //----- Merge data headers 
  MergeZSliceHeaders();


}

//-------------------------------------------------------------
void DicomDetectorConstruction::ReadPatientDataFile(const G4String& fname)
{
#ifdef G4VERBOSE
  G4cout << " DicomDetectorConstruction::ReadPatientDataFile opening file " << fname << G4endl;
#endif 
  std::ifstream fin(fname.c_str(), std::ios_base::in);
  if( !fin.is_open() ) {
    G4Exception("DicomDetectorConstruction::ReadPatientDataFil. File not found " + fname );
  }
  //----- Define density differences (maximum density difference to create a new material)
  G4double densityDiff = 0.1; 
  std::map<G4int,G4double>  fDensityDiffs; // to be able to use a different densityDiff for each material
  for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ii++ ){
    fDensityDiffs[ii] = densityDiff; //currently all materials with same difference
  }

  //----- Read data header
  DicomPatientZSliceHeader* sliceHeader = new DicomPatientZSliceHeader( fin );
  fZSliceHeaders.push_back( sliceHeader );

  //----- Read material indices
  G4int nVoxels = sliceHeader->GetNoVoxels();

  //--- If first slice, initiliaze fMateIDs
  if( fZSliceHeaders.size() == 1 ) {
    //fMateIDs = new unsigned int[fNoFiles*nVoxels];
    fMateIDs = new size_t[fNoFiles*nVoxels];

  }

  unsigned int mateID;
  G4int voxelCopyNo = (fZSliceHeaders.size()-1)*nVoxels; // number of voxels from previously read slices
  for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
    fin >> mateID;
    fMateIDs[voxelCopyNo] = mateID;
  }

  //----- Read material densities and build new materials if two voxels have same material but its density is in a different density interval (size of density intervals defined by densityDiff)
  G4double density;
  voxelCopyNo = (fZSliceHeaders.size()-1)*nVoxels; // number of voxels from previously read slices

  for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ )
  {
		fin >> density;

		//-- Get material from list of original materials

		int mateID = fMateIDs[voxelCopyNo];

		// changed  ML2 start

		//	a search over the materials names is performed to avoide the necessity to respect the material name order.

		//	G4Material* mateOrig  = fOriginalMaterials[mateID];  commented  ML2

		G4Material* mateOrig=0; 
		for (int jMat=0;jMat<(int)fOriginalMaterials.size(); jMat++)
		{
			if (fOriginalMaterials[jMat]->GetName()==this->inputMatNames[mateID])
			{
				mateOrig  = fOriginalMaterials[jMat];
				break;
			}
		}

		// end

		//-- Get density bin: middle point of the bin in which the current density is included 
		G4double densityBin = fDensityDiffs[mateID] * (G4int(density/fDensityDiffs[mateID])+0.5);
		//-- Build the new material name 
		G4String newMateName = mateOrig->GetName()+"__"+ftoa((G4float)densityBin);
	    
		//-- Look if a material with this name is already created (because a previous voxel was already in this density bin)
		unsigned int im;
		for( im = 0; im < fMaterials.size(); im++ )
		{
			if( fMaterials[im]->GetName() == newMateName ) 
			{break;}
	    }
		//-- If material is already created use index of this material
		if( im != fMaterials.size() ) 
		{
			fMateIDs[voxelCopyNo] = im;
			//-- else, create the material 
		} else 
		{
			fMaterials.push_back( BuildMaterialWithChangingDensity( mateOrig, densityBin, newMateName ) );
			fMateIDs[voxelCopyNo] = fMaterials.size()-1;
			this->matNames[this->nMatNames]=newMateName;
			this->nMatNames++;
		}
	}

}


//-------------------------------------------------------------
void DicomDetectorConstruction::MergeZSliceHeaders()
{
  //----- Images must have the same dimension ... 
  fZSliceHeaderMerged = new DicomPatientZSliceHeader( *fZSliceHeaders[0] );
  for( unsigned int ii = 1; ii < fZSliceHeaders.size(); ii++ ) {
    *fZSliceHeaderMerged += *fZSliceHeaders[ii];
  };
}

//-------------------------------------------------------------
G4Material* DicomDetectorConstruction::BuildMaterialWithChangingDensity( const G4Material* origMate, G4double density, G4String newMateName )
{
  //----- Copy original material, but with new density
  G4int nelem = origMate->GetNumberOfElements();
  G4Material* mate = new G4Material( newMateName, density*g/cm3, nelem, kStateUndefined, STP_Temperature );

  for( G4int ii = 0; ii < nelem; ii++ ){
    G4double frac = origMate->GetFractionVector()[ii];
    G4Element* elem = const_cast<G4Element*>(origMate->GetElement(ii));
    mate->AddElement( elem, frac );
  }

  return mate;
}

//-----------------------------------------------------------------------
G4String DicomDetectorConstruction::ftoa(float flo)
{
  char ctmp[100];
  gcvt( flo, 10, ctmp );
  return G4String(ctmp);
}


//-------------------------------------------------------------

void DicomDetectorConstruction::ConstructPatientContainer()
{
  //---- Extract number of voxels and voxel dimensions
  nVoxelX = fZSliceHeaderMerged->GetNoVoxelX();
  nVoxelY = fZSliceHeaderMerged->GetNoVoxelY();
  nVoxelZ = fZSliceHeaderMerged->GetNoVoxelZ();

  voxelHalfDimX = fZSliceHeaderMerged->GetVoxelHalfX();
  voxelHalfDimY = fZSliceHeaderMerged->GetVoxelHalfY();
  voxelHalfDimZ = fZSliceHeaderMerged->GetVoxelHalfZ();
#ifdef G4VERBOSE
  G4cout << " nVoxelX " << nVoxelX << " voxelHalfDimX " << voxelHalfDimX <<G4endl;
  G4cout << " nVoxelY " << nVoxelY << " voxelHalfDimY " << voxelHalfDimY <<G4endl;
  G4cout << " nVoxelZ " << nVoxelZ << " voxelHalfDimZ " << voxelHalfDimZ <<G4endl;
  G4cout << " totalPixels " << nVoxelX*nVoxelY*nVoxelZ <<  G4endl;
#endif

  //----- Define the volume that contains all the voxels
  container_solid = new G4Box("PhantomContainer",nVoxelX*voxelHalfDimX,nVoxelY*voxelHalfDimY,nVoxelZ*voxelHalfDimZ);
  container_logic = 
    new G4LogicalVolume( container_solid, 
			 fMaterials[0],  //the material is not important, it will be fully filled by the voxels
			 "PhantomContainer", 
			 0, 0, 0 );
  //--- Place it on the world
//   G4double offsetX = (fZSliceHeaderMerged->GetMaxX() + fZSliceHeaderMerged->GetMinX() ) /2.;
//   G4double offsetY = (fZSliceHeaderMerged->GetMaxY() + fZSliceHeaderMerged->GetMinY() ) /2.;
//   G4double offsetZ = (fZSliceHeaderMerged->GetMaxZ() + fZSliceHeaderMerged->GetMinZ() ) /2.;
//   G4ThreeVector posCentreVoxels(offsetX,offsetY,offsetZ);

G4ThreeVector posCentreVoxels(0.,0.,0.); // ML2 

#ifdef G4VERBOSE
  G4cout << " placing voxel container volume at " << posCentreVoxels << G4endl;
#endif
		G4RotationMatrix *rotation=new G4RotationMatrix();

		rotation->rotateZ(180.*deg);
  container_phys = 
    new G4PVPlacement(rotation,  // rotation
		      posCentreVoxels,
		      container_logic,     // The logic volume
		      "PhantomContainer",  // Name
		      world_logic,  // Mother
		      false,           // No op. bool.
		      1);              // Copy number


}


#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//-------------------------------------------------------------
void DicomDetectorConstruction::SetScorer(G4LogicalVolume* voxel_logic)
{
	voxel_logic->SetSensitiveDetector(this->sensDet);

	this->regVol->AddRootLogicalVolume(voxel_logic);
	voxel_logic->SetRegion(this->regVol);
}

