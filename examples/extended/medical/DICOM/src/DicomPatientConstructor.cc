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
#include "DicomPatientConstructor.hh"
#include "DicomGeometry.hh"
#include "DicomPatientParameterisation.hh"
#include "DicomConfiguration.hh"
#include "DicomPrimaryGeneratorAction.hh"

#include <fstream>
#include <string.h>
#include <stdio.h>
#include <math.h>

G4int DicomPatientConstructor::FindingNbOfVoxels(G4double MaxDensity , G4double MinDensity)
{
  FILE* lecturepref = fopen("Data.dat","r");
  fscanf(lecturepref,"%s",compressionbuf);
  compression = atoi(compressionbuf);
  fscanf(lecturepref,"%s",maxbuf);
  G4int max = atoi(maxbuf);    
  G4int copy_counter = 0;
  for (G4int i=1;i<=max;i++)
    {
      fscanf(lecturepref,"%s",name);
      sprintf(fullname,"%s.g4",name);
      readData = fopen(fullname,"r");
      fscanf(readData,"%s %s",rowsbuf,columnsbuf);
      G4int rows=atoi(rowsbuf);
      G4int columns=atoi(columnsbuf);
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
      for (G4int j=1;j<=lenr;j++)
        {
	  for (G4int w=1;w<=lenc;w++)
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
// This function reads the Contour information included in Plan.roi
// (this file is a set of region of interest for each dicom images made on ADAC pinnacle treatment planning software)
// This code is still under construction...
void DicomPatientConstructor::readContour()
{

  char ROIplanLine[2000];
  char Word_1[300],Word_2[300],Word_3[300],Word_4[300],Word_5[300];
  FILE* readingContours;

  // the Flag is set for a certain data acquisition
  // 	Flag = 1 -> contours definition, lenROI
  // 	Flag = 2 -> points definition, lenPOINTS
  // 	Flag = 3 -> curves definition, lenCURVE
  G4int Flag=0, lenROI=0, lenPOINTS=0, lenPOINTSref=0, lenCURVE=0, lenCURVEref=0;
  lenCURVE = 0;
  lenROI = 0;
  lenPOINTS = 0;

  readingContours = fopen("Plan.roi","r");
  if ( (G4int *)readingContours == 0 )
    {
      printf("### No contours file ('Plan.roi')\n");
      flag_contours=0;
    }
  else if ( (G4int *)readingContours != 0 )
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
	      G4int x=0;
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
G4bool DicomPatientConstructor::isWithin(G4double X , G4double Y , G4double Z )
{
  G4int x,i,j;
  G4int state=0;  // This variable = 0 if the point is outside curve
  G4int lenPOINTS=-1;  // Number of points avaible in a curve, to be initialised each time!
  for(x=1;x<=MaxCurve;x++)
    {
      lenPOINTS=(G4int)(ContoursX[x][0]);  // The number of points in the curve in the
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

  return (G4bool)state;  // current state false=out, true=inside
}

