//   $tigre.2@sympatico.ca, louis.archambault@phy.ulaval.ca
//   06/12/02

//*******************************************************
//
// DicomPatientConstructor.cc :
//	- Initialisation of the construction of DICM images
//	- Reading contour information included in Plan.roi
//	  (Region of interest) *** NOT FULLY WORKING YET ***
//	- Definitions are in DicomGeometry.hh
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
#include "DicomConfiguration.hh"
#include "DicomPrimaryGeneratorAction.hh"

#include <fstream>
#include <string.h>
#include <stdio.h>
#include <math.h>

using namespace std;

void DicomGeometry::patientConstruction()
{
  DicomConfiguration* ReadConfiguration = new DicomConfiguration;
  ReadConfiguration->ReadDataFile();					// images must have the same dimension
  ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[0] );		//  open a .g4 file to read some values
		
  PatientX = (ReadConfiguration->CompressionUsed*(ReadConfiguration->GetXPixelSpacing())/2.0) *mm;
  PatientY = (ReadConfiguration->CompressionUsed*(ReadConfiguration->Y_PixelSpacing)/2.0) *mm;
  PatientZ = ((ReadConfiguration->SliceTickness/2.0) *mm);

  // Logical Box to place Parameteristion inside it
  Attributes_param = new G4VisAttributes();
  Attributes_param->SetForceSolid(false);
  Attributes_param->SetColour(red=1.,green=0.,blue=0.,alpha=1.);
  G4int totalNumberOfFile = ReadConfiguration->GetTotalNumberOfFile(); 
  G4int totalRows = ReadConfiguration->GetTotalRows(); 
  G4int totalColumns = ReadConfiguration->GetTotalColumns();
  Parameterisation_Box = new G4Box("Parameterisation Mother", totalColumns*(ReadConfiguration->GetXPixelSpacing())/2.*mm, totalRows*(ReadConfiguration->Y_PixelSpacing)/2.*mm, totalNumberOfFile*(ReadConfiguration->SliceTickness)/2.*mm);
  logical_param = new G4LogicalVolume(Parameterisation_Box,air,"Parameterisation Mother (logical)");
  logical_param->SetVisAttributes(Attributes_param);

  double MiddleLocationValue=0;
  for (int i=0;i< totalNumberOfFile;i++)
    {
      ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[i] );
      MiddleLocationValue=MiddleLocationValue+ReadConfiguration->SliceLocation;
    }
  MiddleLocationValue=MiddleLocationValue/totalNumberOfFile;
    
  G4ThreeVector origin( 0.*mm,0.*mm,MiddleLocationValue*mm );
  physical_param =  new G4PVPlacement(0,origin,logical_param,"Parameterisation Mother (placement)",logicWorld,false,0);

  //  Parametrisation of Patient put inside
  LungINhale = new G4Box("LungINhale",PatientX,PatientY,PatientZ);
  Logical_LungINhale = new G4LogicalVolume(LungINhale,lunginhale,"Logical_LungINhale",0,0,0);
  //Logical_LungINhale->SetVisAttributes(Attributes_LungINhale);
  Param_LungINhale = new DicomPatientParameterisation(FindingNbOfVoxels(2.0 , 0.207),  2.0 , 0.207 ,
						      lunginhale,
						      lungexhale,
						      adipose_tissue,
						      breast,
						      phantom,
						      muscle,
						      liver,
						      dense_bone,
						      trabecular_bone);
  Physical_LungINhale = new G4PVParameterised( "Physical_LungINhale" , Logical_LungINhale, logical_param, kZAxis, FindingNbOfVoxels(2.0 , 0.207), Param_LungINhale );
}

int DicomGeometry::FindingNbOfVoxels(double MaxDensity , double MinDensity)
{

  FILE* lecturepref = fopen("Data.dat","r");
  fscanf(lecturepref,"%s",compressionbuf);
  compression=atoi(compressionbuf);
  fscanf(lecturepref,"%s",maxbuf);
  max=atoi(maxbuf);    
  int copy_counter=0;
  for (int i=1;i<=max;i++)
    {
      fscanf(lecturepref,"%s",name);
      sprintf(fullname,"%s.g4",name);
      readData =  fopen(fullname,"r");

      fscanf(readData,"%s %s",rowsbuf,columnsbuf);
      rows=atoi(rowsbuf);
      columns=atoi(columnsbuf);
      fscanf(readData,"%s %s",pixel_spacing_Xbuf,pixel_spacing_Ybuf);
      pixel_spacing_X=atof(pixel_spacing_Xbuf);
      pixel_spacing_Y=atof(pixel_spacing_Ybuf);
      fscanf(readData,"%s",SliceTicknessbuf);
      SliceTickness=atoi(SliceTicknessbuf);
      fscanf(readData,"%s",SliceLocationbuf);
      SliceLocation=atof(SliceLocationbuf);
      fscanf(readData,"%s",compressionbuf);
      compression=atoi(compressionbuf);
      lenr=abs(rows/compression);
      lenc=abs(columns/compression);
      for (int j=1;j<=lenr;j++)
        {
	  for (int w=1;w<=lenc;w++)
            {
	      if ( fscanf(readData,"%s",Densitybuf) != -1 )
                {
		  if ( atof(Densitybuf) >= MinDensity && atof(Densitybuf) <= MaxDensity )
                    {
		      Density.push_back( atof(Densitybuf) );
		      copy_counter++;
                    }
                }
            }
        }
      fclose(readData);
    }
  return copy_counter;
}

void DicomGeometry::InitialisationOfMaterials()
{
  // Creating elements :
  G4double z,a,density;
  G4String name, symbol;
  elC   = new G4Element(name="Carbon",symbol="C",z=6.0,a=12.011 * g/mole);
  elH   = new G4Element(name="Hydrogen",symbol="H",z=1.0,a=1.008  * g/mole);
  elN   = new G4Element(name="Nitrogen",symbol="N",z=7.0,a=14.007 * g/mole);
  elO   = new G4Element(name="Oxygen",symbol="O",z=8.0,a=16.00  * g/mole);
  elNa  = new G4Element(name="Sodium",symbol="Na",z=11.0,a=22.98977* g/mole);
  elS   = new G4Element(name="Sulfur",symbol="S",z=16.0,a=32.065* g/mole);
  elCl  = new G4Element(name="Chlorine",symbol="P",z=17.0,a=35.453* g/mole);
  elK   = new G4Element(name="Potassium",symbol="P",z=19.0,a=30.0983* g/mole);
  elP   = new G4Element(name="Phosphorus",symbol="P",z=30.0,a=30.973976* g/mole);
  elFe  = new G4Element(name="Iron",symbol="Fe",z=26,a=56.845* g/mole);
  elMg  = new G4Element(name="Magnesium",symbol="Mg",z=12.0,a=24.3050* g/mole);
  elCa  = new G4Element(name="Calcium",symbol="Ca",z=20.0,a=40.078* g/mole);

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
  //trabecular_bone->AddElement(elCa,0.074);// not = 1
  trabecular_bone->AddElement(elFe,0.001);

}

// This function reads the Contour information included in Plan.roi
// (this file is a set of region of interest for each dicom images made on ADAC pinnacle treatment planning software)
// This code is still under construction...
void DicomGeometry::readContour()
{

  char ROIplanLine[2000];
  char Word_1[300],Word_2[300],Word_3[300],Word_4[300],Word_5[300];
  FILE* readingContours;

  // the Flag is set for a certain data acquisition
  // 	Flag = 1 -> contours definition, lenROI
  // 	Flag = 2 -> points definition, lenPOINTS
  // 	Flag = 3 -> curves definition, lenCURVE
  int Flag=0, lenROI=0, lenPOINTS=0, lenPOINTSref=0, lenCURVE=0, lenCURVEref=0;
  lenCURVE = 0;
  lenROI = 0;
  lenPOINTS = 0;

  readingContours = fopen("Plan.roi","r");
  if ( (int *)readingContours == 0 )
    {
      printf("### No contours file ('Plan.roi')\n");
      flag_contours=0;
    }
  else if ( (int *)readingContours != 0 )
    {
      flag_contours=1;
      printf("### There is a contour file ('Plan.roi')\n");
      while( fgets(ROIplanLine,2000,readingContours) )
        {
	  if ( ROIplanLine[0] == '/' && ROIplanLine[1] == '/' )
            {
	      // This is a comments line, nothing to do with it
            }
	  else if ( Flag == 1 )
            {
	      sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      // we seek :: num_curve = 39;
	      if( Word_1[0] == 'n' &&  Word_1[1] == 'u' &&  Word_1[2] == 'm' &&  Word_1[3] == '_' &&  Word_1[4] == 'c' &&  Word_1[5] == 'u' &&  Word_1[6] == 'r' &&  Word_1[7] == 'v' &&  Word_1[8] == 'e' )
                {
		  lenCURVEref = atoi(Word_3);
		  Flag=0;
                }
            }
	  else if ( Flag == 2 )
            {
	      sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      // we seek :: num_points = 131
	      if( Word_1[0] == 'n' &&  Word_1[1] == 'u' &&  Word_1[2] == 'm' &&  Word_1[3] == '_' &&  Word_1[4] == 'p' &&  Word_1[5] == 'o' &&  Word_1[6] == 'i' &&  Word_1[7] == 'n' &&  Word_1[8] == 't' )
                {
		  lenPOINTSref = atoi(Word_3);
		  Flag=0;
                }
            }
	  else if ( Flag == 3 )
            {
	      sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      if ( lenPOINTS == 0 )
                {
		  ContoursX[lenCURVE][lenPOINTS]=lenPOINTSref;
		  ContoursY[lenCURVE][lenPOINTS]=lenPOINTSref;
		  ContoursZ[lenCURVE][lenPOINTS]=lenPOINTSref;
                }
	      lenPOINTS++;
	      if( Word_1[0] == '}' &&  Word_1[1] == ';')
		Flag=0;
	      else
                {
		  ContoursX[lenCURVE][lenPOINTS]=atof(Word_1);
		  ContoursY[lenCURVE][lenPOINTS]=atof(Word_2);
		  ContoursZ[lenCURVE][lenPOINTS]=atof(Word_3);
                }
            }
	  else
            {
	      int x=0;
	      sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      while( Word_1[x] != 0 )
                {
		  x++;
		  if ( Word_1[x] == '=' && Word_1[x+1] == '{')
                    {
		      if ( Word_1[x-3] == 'r' &&  Word_1[x-2] == 'o' && Word_1[x-1] == 'i' )
                        {
			  Flag = 1;
			  lenROI++;
                        }
		      else if ( Word_1[x-5] == 'c' &&  Word_1[x-4] == 'u' && Word_1[x-3] == 'r' && Word_1[x-2] == 'v' && Word_1[x-1] == 'e' )
                        {
			  Flag = 2;
			  lenCURVE++;
			  MaxCurve=lenCURVE;
                        }
		      else if ( Word_1[x-6] == 'p' &&  Word_1[x-5] == 'o' && Word_1[x-4] == 'i' && Word_1[x-3] == 'n' && Word_1[x-2] == 't' && Word_1[x-1] == 's' )
                        {
			  Flag = 3;
			  lenPOINTS=0;
                        }
		      else if ( x == 2000 )
                        {
			  printf("### ERROR, No definition found !\n");
			  printf("%s\n",Word_1);
                        }
                    }
                }
            }
        }

      if ( (int *)readingContours != 0 )
        {
	  fclose(readingContours);
        }
    }
}

// This function tell if the coordinates (X,Y,Z) is within one of the curves from plan.roi
bool DicomGeometry::isWithin(double X , double Y , double Z )
{
  int x,i,j;
  int state=0;  // This variable = 0 if the point is outside curve
  int lenPOINTS=-1;  // Number of points avaible in a curve, to be initialised each time!
  for(x=1;x<=MaxCurve;x++)
    {
      lenPOINTS=(int)(ContoursX[x][0]);  // The number of points in the curve in the
      // first element of the array
      if (  ContoursZ[x][1] ==  Z )  // Z[x][1] is the same for all Z
        {
	  for (i = 1, j = lenPOINTS; i <= lenPOINTS; j = i++)
            {
	      if (	( ( (ContoursY[x][i]<=Y) && (Y<ContoursY[x][j]) ) || ( (ContoursY[x][j]<=Y) && (Y<ContoursY[x][i]) ) )	&&
                        (X < (ContoursX[x][j] - ContoursX[x][i])*(Y - ContoursY[x][i]) / (ContoursY[x][j] - ContoursY[x][i]) + ContoursX[x][i])	)  // If this is true we've just pass a limit of the curve
                {
		  state=!state;  // the state changes
                }
            }
        }

    }

  return (bool)state;  // current state false=out, true=inside
}

