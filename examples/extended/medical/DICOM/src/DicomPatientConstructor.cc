//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************

//*******************************************************
//
// DicomPatientConstructor.cc :
//	- Initialisation of the construction of DICM images
//	- Reading contour information included in Plan.roi
//	  (Region of interest) *** NOT FULLY WORKING YET ***
//	- Definitions are in DicomGeometry.hh
//
//*******************************************************

#include "DicomPatientConstructor.hh"

#include "globals.hh"
#include "G4String.hh"

G4int DicomPatientConstructor::FindingNbOfVoxels(G4double maxDensity , G4double minDensity)
{
  FILE* lecturePref = std::fopen("Data.dat","r");
  FILE* readData;
  char  maxBuffer[300];
  char compressionBuffer[300];
  char name[300];
  std::fscanf(lecturePref,"%s",compressionBuffer);
  G4int compression;
  compression = atoi(compressionBuffer);

  std::fscanf(lecturePref,"%s",maxBuffer);
  G4int max = atoi(maxBuffer);    
  G4int copyCount = 0;

  for ( G4int i = 1;i <= max;i++ )
    {
      char fullname[300];
      std::fscanf(lecturePref,"%s",name);
      std::sprintf(fullname,"%s.g4",name);
      readData = std::fopen(fullname,"r");
 
      char rowsBuffer[300],columnsBuffer[300];
      std::fscanf(readData,"%s %s",rowsBuffer,columnsBuffer);
      G4int rows = atoi(rowsBuffer);
      G4int columns = atoi(columnsBuffer);

      char pixelSpacingXBuffer[300],pixelSpacingYBuffer[300];
      std::fscanf(readData,"%s %s",pixelSpacingXBuffer,pixelSpacingYBuffer);
      pixelSpacingX = atof(pixelSpacingXBuffer);
      pixelSpacingY = atof(pixelSpacingYBuffer);
      
      char sliceThicknessBuf[300];
      std::fscanf(readData,"%s",sliceThicknessBuf);
      sliceThickness = atoi(sliceThicknessBuf);

      char sliceLocationBuf[300];
      std::fscanf(readData,"%s",sliceLocationBuf);
      sliceLocation = atof(sliceLocationBuf);

      std::fscanf(readData,"%s",compressionBuffer);
      compression = atoi(compressionBuffer);
      lenr = abs(rows/compression);
      lenc = abs(columns/compression);
      char densityBuffer[300];
      std::vector<G4double> density;
      for ( G4int j = 1;j <= lenr;j++ )
        {
	  for ( G4int w = 1;w <= lenc;w++ )
            {
	      if ( std::fscanf(readData,"%s",densityBuffer) != -1 )
                {
		  if ( atof(densityBuffer) >= minDensity && atof(densityBuffer) <= maxDensity )
                    {
		      density.push_back( atof(densityBuffer) );
		      copyCount++;
                    }
                }
            }
        }
      std::fclose(readData);
    }
  return copyCount;
}
 
