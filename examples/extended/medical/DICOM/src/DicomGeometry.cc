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

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "DicomGeometry.hh"
#include "DicomPatientParameterisation.hh"
#include "DicomPatientConstructor.hh"
#include "DicomConfiguration.hh"


DicomGeometry::DicomGeometry()
{
  patientConstructor = new DicomPatientConstructor();
 
  solidWorld = 0;
  logicWorld = 0;
  physiWorld = 0;
  parameterisedPhysVolume = 0;
  physicalLungINhale = 0;
}

DicomGeometry::~DicomGeometry()
{ }

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
  G4Material* trabecularBone = new G4Material( "Skeleton_Spongiosa", 
					       density = 1159.*kg/m3, 
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
  G4Material* denseBone = new G4Material( "Skeleton_ribs", 
					  density = 1575.*kg/m3, 
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
  G4Material* liver = new G4Material( "Liver", 
				      density = 1071.*kg/m3, 
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
  G4Material* muscle = new G4Material( "Muscle", 
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
  G4Material* phantom = new G4Material( "Phantom", 
					density = 1.018*kg/m3, 
					numberofElements = 2 );
  phantom->AddElement(elH,0.112);
  phantom->AddElement(elO,0.888);     


  // Breast
  G4Material* breast = new G4Material( "Breast", 
                           density = 990.*kg/m3, 
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
  G4Material* adiposeTissue = new G4Material( "adipose_tissue", 
                                  density = 967.*kg/m3, 
                                  numberofElements = 7);
  adiposeTissue->AddElement(elH,0.114);
  adiposeTissue->AddElement(elC,0.598);
  adiposeTissue->AddElement(elN,0.007);
  adiposeTissue->AddElement(elO,0.278);
  adiposeTissue->AddElement(elNa,0.001);
  adiposeTissue->AddElement(elS,0.001);
  adiposeTissue->AddElement(elCl,0.001);


  //  LungEXhale
  G4Material* lungExhale = new G4Material( "Lung_Exhale", 
					   density = 508.*kg/m3, 
					   numberofElements = 9 );
  lungExhale->AddElement(elH,0.103);
  lungExhale->AddElement(elC,0.105);
  lungExhale->AddElement(elN,0.031);
  lungExhale->AddElement(elO,0.749);
  lungExhale->AddElement(elNa,0.002);
  lungExhale->AddElement(elP,0.002);
  lungExhale->AddElement(elS,0.003);
  lungExhale->AddElement(elCl,0.002);
  lungExhale->AddElement(elK,0.003);

  //  LungINhale
  G4Material* lungInhale = new G4Material( "Lung_Inhale", 
                               density = 217.*kg/m3, 
                               numberofElements = 9);
  lungInhale->AddElement(elH,0.103);
  lungInhale->AddElement(elC,0.105);
  lungInhale->AddElement(elN,0.031);
  lungInhale->AddElement(elO,0.749);
  lungInhale->AddElement(elNa,0.002);
  lungInhale->AddElement(elP,0.002);
  lungInhale->AddElement(elS,0.003);
  lungInhale->AddElement(elCl,0.002);
  lungInhale->AddElement(elK,0.003);


  // Air
  G4Material* air = new G4Material( "Air",
				    1.290*mg/cm3,
				    numberofElements = 2 );
  air->AddElement(elN, 0.7);
  air->AddElement(elO, 0.3); 
}

void DicomGeometry::PatientConstruction()
{
  DicomConfiguration readConfiguration;
  readConfiguration.ReadDataFile();
		
  G4String listOfFile = readConfiguration.GetListOfFile()[0];	
  // images must have the same dimension ... 
  readConfiguration.ReadG4File( listOfFile );
  // open a .g4 file to read some values ...
  
  G4int compressionUsed = readConfiguration.IsCompressionUsed();	
  G4double sliceThickness = readConfiguration.GetSliceThickness();
  G4double xPixelSpacing = readConfiguration.GetXPixelSpacing(); 
  G4double yPixelSpacing = readConfiguration.GetYPixelSpacing(); 	
  G4int totalNumberOfFile = readConfiguration.GetTotalNumberOfFile(); 
  G4int totalRows = readConfiguration.GetTotalRows(); 
  G4int totalColumns = readConfiguration.GetTotalColumns();

  G4double patientX = (compressionUsed*(xPixelSpacing)/2.0) *mm;
  G4double patientY = (compressionUsed*(yPixelSpacing)/2.0) *mm;
  G4double patientZ = ((sliceThickness/2.0) *mm);

  G4VisAttributes* visualisationAttribute = new G4VisAttributes();
  visualisationAttribute->SetForceSolid(false);
  visualisationAttribute->SetColour( 1., 0., 0., 1. );
  
  //Building up the parameterisation ...
  G4Box* parameterisedBox = new G4Box( "Parameterisation Mother", 
				       totalColumns * (xPixelSpacing ) / 2.* mm, 
				       totalRows * (yPixelSpacing) / 2.* mm,
				       totalNumberOfFile * (sliceThickness) / 2.* mm);
  G4Material* air = G4Material::GetMaterial("Air");
  G4LogicalVolume* parameterisedLogicalvolume = 
    new G4LogicalVolume( parameterisedBox,
			 air,
			 "Parameterisation Mother (logical)" );
  parameterisedLogicalvolume->SetVisAttributes(visualisationAttribute);
  
  G4double middleLocationValue = 0;
  for ( G4int i=0; i< totalNumberOfFile;i++ )
    {
      readConfiguration.ReadG4File( readConfiguration.GetListOfFile()[i] );
      middleLocationValue = middleLocationValue + readConfiguration.GetSliceLocation();
    }
  middleLocationValue = middleLocationValue / totalNumberOfFile;
    
  G4ThreeVector origin( 0.*mm, 0.*mm, middleLocationValue * mm );
    parameterisedPhysVolume =  new G4PVPlacement( 0,
                                                  origin,
		                                  parameterisedLogicalvolume,
		                                  "Parameterisation Mother (placement)",
		                                  logicWorld,
		                                  false,
		                                  0 );

  G4Box* LungINhale = new G4Box( "LungINhale", patientX, patientY, patientZ);

  G4Material* lungInhaleMaterial = G4Material::GetMaterial("Lung_Inhale");
  G4LogicalVolume* logicLungInHale = new G4LogicalVolume(LungINhale,lungInhaleMaterial,"Logical_LungINhale",0,0,0);

  // ---- MGP ---- Numbers (2.0, 0.207) to be removed from code; move to const
  G4int numberOfVoxels = patientConstructor->FindingNbOfVoxels(2.0,0.207);

  G4VPVParameterisation* paramLungINhale = new DicomPatientParameterisation
    ( numberOfVoxels,
      2.0 , 0.207 );

  physicalLungINhale = 
    new G4PVParameterised( "Physical_LungINhale" , 
			   logicLungInHale, 
			   parameterisedLogicalvolume,
			   kZAxis, numberOfVoxels, 
			   paramLungINhale );
}

G4VPhysicalVolume* DicomGeometry::Construct()
{
  InitialisationOfMaterials();

  G4double worldXDimension = 1.*m;
  G4double worldYDimension = 1.*m;
  G4double worldZDimension = 1.*m;

  solidWorld = new G4Box( "WorldSolid",
                          worldXDimension,
                          worldYDimension,
                          worldZDimension );

  G4Material* air = G4Material::GetMaterial("Air");

  logicWorld = new G4LogicalVolume( solidWorld, 
                                    air, 
                                    "WorldLogical", 
                                    0, 0, 0 );

  physiWorld = new G4PVPlacement( 0,
                                  G4ThreeVector(0,0,0),
                                  "World",
                                  logicWorld,
                                  0,
                                  false,
                                  0 );
  PatientConstruction();

  return physiWorld;
}
