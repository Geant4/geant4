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

#include "DicomPatientConstructor.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4String.hh"
#include <stdio.h>
#include <math.h>

G4int DicomPatientConstructor::FindingNbOfVoxels(G4double maxDensity , G4double minDensity)
{
  FILE* lecturePref = G4std::fopen("Data.dat","r");

  G4std::fscanf(lecturePref,"%s",compressionBuf);
  compression = atoi(compressionBuf);

  G4std::fscanf(lecturePref,"%s",maxBuf);
  G4int max = atoi(maxBuf);    
  G4int copy_counter = 0;

  for (G4int i=1;i<=max;i++)
    {
      G4std::fscanf(lecturePref,"%s",name);
      G4std::sprintf(fullname,"%s.g4",name);
      readData = G4std::fopen(fullname,"r");
 
      G4std::fscanf(readData,"%s %s",rowsBuf,columnsBuf);
      G4int rows = atoi(rowsBuf);
      G4int columns = atoi(columnsBuf);

      G4std::fscanf(readData,"%s %s",pixelSpacingXBuf,pixelSpacingYBuf);
      pixelSpacingX = atof(pixelSpacingXBuf);
      pixelSpacingY = atof(pixelSpacingYBuf);

      G4std::fscanf(readData,"%s",sliceTicknessBuf);
      sliceTickness = atoi(sliceTicknessBuf);

      G4std::fscanf(readData,"%s",sliceLocationBuf);
      sliceLocation = atof(sliceLocationBuf);

      G4std::fscanf(readData,"%s",compressionBuf);
      compression = atoi(compressionBuf);
      lenr = abs(rows/compression);
      lenc = abs(columns/compression);
 
      for (G4int j=1;j<=lenr;j++)
        {
	  for (G4int w=1;w<=lenc;w++)
            {
	      if ( G4std::fscanf(readData,"%s",densityBuf) != -1 )
                {
		  if ( atof(densityBuf) >= minDensity && atof(densityBuf) <= maxDensity )
                    {
		      density.push_back( atof(densityBuf) );
		      copy_counter++;
                    }
                }
            }
        }
      G4std::fclose(readData);
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
  // 	flag = 1 -> contours definition, lenROI
  // 	flag = 2 -> points definition, lenPoints
  // 	flag = 3 -> curves definition, lenCurve
  G4int flag=0, lenROI=0, lenPoints=0, lenPointsRef=0, lenCurve=0, lenCurveRef=0;
  lenCurve = 0;
  lenROI = 0;
  lenPoints = 0;

  readingContours = G4std::fopen("Plan.roi","r");
  if ( (G4int *)readingContours == 0 )
    {
      G4std::printf("### No contours file ('Plan.roi')\n");
      flag_contours = 0;
    }
  else if ( (G4int *)readingContours != 0 )
    {
      flag_contours = 1;
      G4std::printf("### There is a contour file ('Plan.roi')\n");
      while( G4std::fgets(ROIplanLine,2000,readingContours) )
        {
	  if ( ROIplanLine[0] == '/' && ROIplanLine[1] == '/' )
            {
	      // This is a comments line, nothing to do with it
            }
	  else if ( flag == 1 )
            {
	      G4std::sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      // we seek :: num_curve = 39;
	      if( Word_1[0] == 'n' &&  Word_1[1] == 'u' &&  Word_1[2] == 'm' &&  Word_1[3] == '_' &&  Word_1[4] == 'c' &&  Word_1[5] == 'u' &&  Word_1[6] == 'r' &&  Word_1[7] == 'v' &&  Word_1[8] == 'e' )
                {
		  lenCurveRef = atoi(Word_3);
		  flag=0;
                }
            }
	  else if ( flag == 2 )
            {
	      G4std::sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      // we seek :: num_points = 131
	      if( Word_1[0] == 'n' &&  Word_1[1] == 'u' &&  Word_1[2] == 'm' &&  Word_1[3] == '_' &&  Word_1[4] == 'p' &&  Word_1[5] == 'o' &&  Word_1[6] == 'i' &&  Word_1[7] == 'n' &&  Word_1[8] == 't' )
                {
		  lenPointsRef = atoi(Word_3);
		  flag = 0;
                }
            }
	  else if ( flag == 3 )
            {
	      G4std::sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      if ( lenPoints == 0 )
                {
		  contoursX[lenCurve][lenPoints] = lenPointsRef;
		  contoursY[lenCurve][lenPoints] = lenPointsRef;
		  contoursZ[lenCurve][lenPoints] = lenPointsRef;
                }
	      lenPoints++;
	      if( Word_1[0] == '}' &&  Word_1[1] == ';')
		flag = 0;
	      else
                {
		  contoursX[lenCurve][lenPoints] = atof(Word_1);
		  contoursY[lenCurve][lenPoints] = atof(Word_2);
		  contoursZ[lenCurve][lenPoints] = atof(Word_3);
                }
            }
	  else
            {
	      G4int x = 0;
	      G4std::sscanf(ROIplanLine,"%s %s %s %s %s",Word_1,Word_2,Word_3,Word_4,Word_5);
	      while( Word_1[x] != 0 )
                {
		  x++;
		  if ( Word_1[x] == '=' && Word_1[x+1] == '{')
                    {
		      if ( Word_1[x-3] == 'r' &&  Word_1[x-2] == 'o' && Word_1[x-1] == 'i' )
                        {
			  flag = 1;
			  lenROI++;
                        }
		      else if ( Word_1[x-5] == 'c' &&  Word_1[x-4] == 'u' && Word_1[x-3] == 'r' && Word_1[x-2] == 'v' && Word_1[x-1] == 'e' )
                        {
			  flag = 2;
			  lenCurve++;
			  maxCurve = lenCurve;
                        }
		      else if ( Word_1[x-6] == 'p' &&  Word_1[x-5] == 'o' && Word_1[x-4] == 'i' && Word_1[x-3] == 'n' && Word_1[x-2] == 't' && Word_1[x-1] == 's' )
                        {
			  flag = 3;
			  lenPoints = 0;
                        }
		      else if ( x == 2000 )
                        {
			  G4std::printf("### ERROR, No definition found !\n");
			  G4std::printf("%s\n",Word_1);
                        }
                    }
                }
            }
        }

      if ( (G4int *)readingContours != 0 )
        {
	  G4std::fclose(readingContours);
        }
    }
}

// This function tell if the coordinates (X,Y,Z) is within one of the curves from plan.roi
G4bool DicomPatientConstructor::isWithin(G4double X , G4double Y , G4double Z )
{
  G4int x,i,j;
  G4int state = 0;  // This variable = 0 if the point is outside curve
  G4int lenPoints = -1;  // Number of points avaible in a curve, to be initialised each time!
  for(x=1;x<=maxCurve;x++)
    {
      lenPoints = (G4int)(contoursX[x][0]);  // The number of points in the curve in the
      // first element of the array
      if (  contoursZ[x][1] ==  Z )  // Z[x][1] is the same for all Z
        {
	  for (i = 1, j = lenPoints; i <= lenPoints; j = i++)
            {
	      if (	( ( (contoursY[x][i]<=Y) && (Y<contoursY[x][j]) ) || ( (contoursY[x][j]<=Y) && (Y<contoursY[x][i]) ) )	&&
                        (X < (contoursX[x][j] - contoursX[x][i])*(Y - contoursY[x][i]) / (contoursY[x][j] - contoursY[x][i]) + contoursX[x][i])	)  // If this is true we've just pass a limit of the curve
                {
		  state=!state;  // the state changes
                }
            }
        }

    }

  return (G4bool)state;  // current state false=out, true=inside
}

