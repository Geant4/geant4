//*******************************************************
//
// DicomGeometry.cc :
//	- Starting the building of the geometry
// 	- Creation of the world and other mother volume
//	- Initialisation of patient geometry
// 	- Definitions are in DicomGeometry.hh
//
// The code was written by :
//	Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
// For more information contact :
//	Louis Archambault louis.archambault@phy.ulaval.ca
// or
//	Luc Beaulieu beaulieu@phy.ulaval.ca
//
// Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//*******************************************************
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4PVParameterised.hh"
#include "DicomGeometry.hh"
#include "DicomPatientParameterisation.hh"
#include "DicomPatientConstructor.hh"
#include "DicomConfiguration.hh"
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <math.h>

using namespace std;

DicomGeometry *DicomGeometry::theDetector=0;

DicomGeometry::DicomGeometry()
{
  patientConstructor = new DicomPatientConstructor();
  theDetector=this;
}

DicomGeometry::~DicomGeometry()
{
  theDetector=0;
  delete patientConstructor;
  delete air;
  delete lunginhale;
  delete lungexhale;
  delete adipose_tissue;
  delete breast;
  delete phantom;
  delete muscle;
  delete liver;
  delete dense_bone;
  delete trabecular_bone;
}
void DicomGeometry::PatientConstruction()
{
  DicomConfiguration* ReadConfiguration = new DicomConfiguration;
  ReadConfiguration->ReadDataFile();					// images must have the same dimension
  ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[0] );		//  open a .g4 file to read some values
		
  PatientX = (ReadConfiguration->IsCompressionUsed()*(ReadConfiguration->GetXPixelSpacing())/2.0) *mm;
  PatientY = (ReadConfiguration->IsCompressionUsed()*(ReadConfiguration->GetYPixelSpacing())/2.0) *mm;
  PatientZ = ((ReadConfiguration->GetSliceThickness()/2.0) *mm);

  // Logical Box to place Parameteristion inside it
  Attributes_param = new G4VisAttributes();
  Attributes_param->SetForceSolid(false);
  Attributes_param->SetColour(red=1.,green=0.,blue=0.,alpha=1.);
  G4int totalNumberOfFile = ReadConfiguration->GetTotalNumberOfFile(); 
  G4int totalRows = ReadConfiguration->GetTotalRows(); 
  G4int totalColumns = ReadConfiguration->GetTotalColumns();
  Parameterisation_Box = new G4Box("Parameterisation Mother", totalColumns*(ReadConfiguration->GetXPixelSpacing())/2.*mm, totalRows*(ReadConfiguration->GetYPixelSpacing())/2.*mm, totalNumberOfFile*(ReadConfiguration->GetSliceThickness())/2.*mm);
  logical_param = new G4LogicalVolume(Parameterisation_Box,air,"Parameterisation Mother (logical)");
  logical_param->SetVisAttributes(Attributes_param);

  G4double MiddleLocationValue=0;
  for (int i=0;i< totalNumberOfFile;i++)
    {
      ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[i] );
      MiddleLocationValue=MiddleLocationValue+ReadConfiguration->GetSliceLocation();
    }
  MiddleLocationValue=MiddleLocationValue/totalNumberOfFile;
    
  G4ThreeVector origin( 0.*mm,0.*mm,MiddleLocationValue*mm );
  physical_param =  new G4PVPlacement(0,origin,logical_param,"Parameterisation Mother (placement)",logicWorld,false,0);

  //  Parametrisation of Patient put inside
  LungINhale = new G4Box("LungINhale",PatientX,PatientY,PatientZ);
  Logical_LungINhale = new G4LogicalVolume(LungINhale,lunginhale,"Logical_LungINhale",0,0,0);
  //Logical_LungINhale->SetVisAttributes(Attributes_LungINhale);
  G4int numberOfVoxels = patientConstructor->FindingNbOfVoxels(2.0,0.207);
  Param_LungINhale = new DicomPatientParameterisation(numberOfVoxels,
                                                       2.0 , 0.207 ,
						      lunginhale,
						      lungexhale,
						      adipose_tissue,
						      breast,
						      phantom,
						      muscle,
						      liver,
						      dense_bone,
						      trabecular_bone);
  Physical_LungINhale = new G4PVParameterised( "Physical_LungINhale" , Logical_LungINhale, logical_param, kZAxis, numberOfVoxels, Param_LungINhale );
}

G4VPhysicalVolume* DicomGeometry::Construct()
{
  InitialisationOfMaterials();

  world_x_width = 1.*m;
  world_y_width = 1.*m;
  world_z_width = 1.*m;

  theWorldDim = G4ThreeVector(world_x_width,world_y_width,world_z_width);

  solidWorld = new G4Box("WorldSolid",world_x_width,world_y_width,world_z_width);
  logicWorld = new G4LogicalVolume( solidWorld, air, "WorldLogical", 0, 0, 0);
  physiWorld = new G4PVPlacement(0,G4ThreeVector(0,0,0),"World",logicWorld,0,false,0);

  patientConstructor -> readContour();  // Contours are not mandatory and are NOT finish yet
  PatientConstruction();

  return physiWorld;
}

void DicomGeometry::InitialisationOfMaterials()
{
  // Creating elements :
  G4double z,a,density;
  G4String name, symbol;

  G4Element* elC = new G4Element(name="Carbon",symbol="C",z=6.0,a=12.011 * g/mole);
  G4Element* elH = new G4Element(name="Hydrogen",symbol="H",z=1.0,a=1.008  * g/mole);
  G4Element* elN = new G4Element(name="Nitrogen",symbol="N",z=7.0,a=14.007 * g/mole);
  G4Element* elO = new G4Element(name="Oxygen",symbol="O",z=8.0,a=16.00  * g/mole);
  G4Element* elNa = new G4Element(name="Sodium",symbol="Na",z=11.0,a=22.98977* g/mole);
  G4Element* elS = new G4Element(name="Sulfur",symbol="S",z=16.0,a=32.065* g/mole);
  G4Element* elCl = new G4Element(name="Chlorine",symbol="P",z=17.0,a=35.453* g/mole);
  G4Element* elK = new G4Element(name="Potassium",symbol="P",z=19.0,a=30.0983* g/mole);
  G4Element* elP = new G4Element(name="Phosphorus",symbol="P",z=30.0,a=30.973976* g/mole);
  G4Element* elFe = new G4Element(name="Iron",symbol="Fe",z=26,a=56.845* g/mole);
  G4Element* elMg = new G4Element(name="Magnesium",symbol="Mg",z=12.0,a=24.3050* g/mole);
  G4Element* elCa = new G4Element(name="Calcium",symbol="Ca",z=20.0,a=40.078* g/mole);

  // Creating Materials :
  G4int nel;

  // Air
  air = new G4Material("Air",1.290*mg/cm3,nel=2);
  air->AddElement(elN, 0.7);
  air->AddElement(elO, 0.3);

  //  LungINhale
  lunginhale = new G4Material("Lung_Inhale", density = 217*kg/m3, nel=9);
  lunginhale->AddElement(elH,0.103);
  lunginhale->AddElement(elC,0.105);
  lunginhale->AddElement(elN,0.031);
  lunginhale->AddElement(elO,0.749);
  lunginhale->AddElement(elNa,0.002);
  lunginhale->AddElement(elP,0.002);
  lunginhale->AddElement(elS,0.003);
  lunginhale->AddElement(elCl,0.002);
  lunginhale->AddElement(elK,0.003);//=1

  //  LungEXhale  (seulement la densite change)
  lungexhale = new G4Material("Lung_Exhale", density = 508*kg/m3, nel=9);
  lungexhale->AddElement(elH,0.103);
  lungexhale->AddElement(elC,0.105);
  lungexhale->AddElement(elN,0.031);
  lungexhale->AddElement(elO,0.749);
  lungexhale->AddElement(elNa,0.002);
  lungexhale->AddElement(elP,0.002);
  lungexhale->AddElement(elS,0.003);
  lungexhale->AddElement(elCl,0.002);
  lungexhale->AddElement(elK,0.003);//=1

  // Adipose tissue
  adipose_tissue = new G4Material("Adipose_tissue", density = 967*kg/m3, nel=7);
  adipose_tissue->AddElement(elH,0.114);
  adipose_tissue->AddElement(elC,0.598);
  adipose_tissue->AddElement(elN,0.007);
  adipose_tissue->AddElement(elO,0.278);
  adipose_tissue->AddElement(elNa,0.001);
  adipose_tissue->AddElement(elS,0.001);
  adipose_tissue->AddElement(elCl,0.001);//=1

  // Breast
  breast = new G4Material("Breast", density = 990*kg/m3, nel=8);
  breast->AddElement(elH,0.109);
  breast->AddElement(elC,0.506);
  breast->AddElement(elN,0.023);
  breast->AddElement(elO,0.358);
  breast->AddElement(elNa,0.001);
  breast->AddElement(elP,0.001);
  breast->AddElement(elS,0.001);
  breast->AddElement(elCl,0.001);	//=1

  // Phantom
  phantom = new G4Material("Phantom", density = 1.018*kg/m3, nel=2);
  phantom->AddElement(elH,0.112);
  phantom->AddElement(elO,0.888);	//=1

  // Muscle
  muscle = new G4Material("Muscle", density = 1061*kg/m3, nel=9);
  muscle->AddElement(elH,0.102);
  muscle->AddElement(elC,0.143);
  muscle->AddElement(elN,0.034);
  muscle->AddElement(elO,0.710);
  muscle->AddElement(elNa,0.001);
  muscle->AddElement(elP,0.002);
  muscle->AddElement(elS,0.003);
  muscle->AddElement(elCl,0.001);
  muscle->AddElement(elK,0.004);	//=1

  // Liver
  liver = new G4Material("Liver", density = 1071*kg/m3, nel=9);
  liver->AddElement(elH,0.102);
  liver->AddElement(elC,0.139);
  liver->AddElement(elN,0.030);
  liver->AddElement(elO,0.716);
  liver->AddElement(elNa,0.002);
  liver->AddElement(elP,0.003);
  liver->AddElement(elS,0.003);
  liver->AddElement(elCl,0.002);
  liver->AddElement(elK,0.003);	//=1

  // Dense Bone (not sure about composition)
  dense_bone = new G4Material("Skeleton_ribs", density = 1575*kg/m3, nel=11);
  dense_bone->AddElement(elH,0.056);
  dense_bone->AddElement(elC,0.235);
  dense_bone->AddElement(elN,0.050);
  dense_bone->AddElement(elO,0.434);
  dense_bone->AddElement(elNa,0.001);
  dense_bone->AddElement(elMg,0.001);
  dense_bone->AddElement(elP,0.072);
  dense_bone->AddElement(elS,0.003);
  dense_bone->AddElement(elCl,0.001);
  dense_bone->AddElement(elK,0.001);
  dense_bone->AddElement(elCa,0.146);
  //dense_bone->AddElement(elCa,0.156);// not = 1

  // Trabecular Bone (not sure about composition)
  trabecular_bone = new G4Material("Skeleton_Spongiosa", density = 1159*kg/m3, nel=12);
  trabecular_bone->AddElement(elH,0.085);
  trabecular_bone->AddElement(elC,0.404);
  trabecular_bone->AddElement(elN,0.058);
  trabecular_bone->AddElement(elO,0.367);
  trabecular_bone->AddElement(elNa,0.001);
  trabecular_bone->AddElement(elMg,0.001);
  trabecular_bone->AddElement(elP,0.034);
  trabecular_bone->AddElement(elS,0.002);
  trabecular_bone->AddElement(elCl,0.002);
  trabecular_bone->AddElement(elK,0.001);
  trabecular_bone->AddElement(elCa,0.044);
  trabecular_bone->AddElement(elFe,0.001);
}



