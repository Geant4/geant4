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

DicomGeometry::DicomGeometry()
{
  patientConstructor = new DicomPatientConstructor();
  trabecularBone = 0;
  denseBone = 0;
  liver = 0;  
  muscle = 0;
  phantom = 0;
  breast = 0;
  adiposeTissue = 0;
  lungexhale = 0;
  lunginhale = 0;
  air = 0;
}

DicomGeometry::~DicomGeometry()
{
  delete air;
  delete lunginhale;
  delete lungexhale;
  delete adiposeTissue;
  delete breast;
  delete phantom;
  delete muscle;
  delete liver;
  delete denseBone;
  delete trabecularBone; 
  delete patientConstructor;
}


void DicomGeometry::InitialisationOfMaterials()
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

  // Trabecular Bone 
  trabecularBone = new G4Material( "Skeleton_Spongiosa", 
                                  density = 1159*kg/m3, 
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
  
  // dense Bone
  denseBone = new G4Material( "Skeleton_ribs", 
                              density = 1575*kg/m3, 
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
 
  // Liver
  liver = new G4Material( "Liver", 
                          density = 1071*kg/m3, 
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
 
  // Muscle
  muscle = new G4Material( "Muscle", 
                           density = 1061*kg/m3, 
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
  
  // Phantom
  phantom = new G4Material( "Phantom", 
                            density = 1.018*kg/m3, 
                            numberofElements = 2 );
  phantom->AddElement(elH,0.112);
  phantom->AddElement(elO,0.888);     


  // Breast
  breast = new G4Material( "Breast", 
                           density = 990*kg/m3, 
                           numberofElements = 8 );
  breast->AddElement(elH,0.109);
  breast->AddElement(elC,0.506);
  breast->AddElement(elN,0.023);
  breast->AddElement(elO,0.358);
  breast->AddElement(elNa,0.001);
  breast->AddElement(elP,0.001);
  breast->AddElement(elS,0.001);
  breast->AddElement(elCl,0.001); 

 // Adipose tissue
  adiposeTissue = new G4Material( "adipose_tissue", 
                                  density = 967*kg/m3, 
                                  numberofElements = 7);
  adiposeTissue->AddElement(elH,0.114);
  adiposeTissue->AddElement(elC,0.598);
  adiposeTissue->AddElement(elN,0.007);
  adiposeTissue->AddElement(elO,0.278);
  adiposeTissue->AddElement(elNa,0.001);
  adiposeTissue->AddElement(elS,0.001);
  adiposeTissue->AddElement(elCl,0.001);

  lungexhale = new G4Material( "Lung_Exhale", 
                               density = 508*kg/m3, 
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

  //  LungINhale
  lunginhale = new G4Material( "Lung_Inhale", 
                               density = 217*kg/m3, 
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


 // Air
  air = new G4Material( "Air",
                        1.290*mg/cm3,
                        numberofElements = 2 );
  air->AddElement(elN, 0.7);
  air->AddElement(elO, 0.3); 
}

void DicomGeometry::PatientConstruction()
{
  DicomConfiguration* ReadConfiguration = new DicomConfiguration;
  ReadConfiguration->ReadDataFile();
					
  // images must have the same dimension ... 
  ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[0] );
  // open a .g4 file to read some values ...
		
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
  for (G4int i=0;i< totalNumberOfFile;i++)
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
						      adiposeTissue,
						      breast,
						      phantom,
						      muscle,
						      liver,
						      denseBone,
						      trabecularBone);
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
