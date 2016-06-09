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
//
//
// DicomConfiguration.cc :
//	- Handling of the header of *.g4 files
// 	- Reading <Data.dat> file

#include "globals.hh"
#include "DicomConfiguration.hh"
#include <fstream>
#include <vector>

G4bool DicomConfiguration::ReadDataFile()
{
  std::ifstream dataFile( "Data.dat" );
  G4String nameOfFileBuffer;
  if ( dataFile.good() != 1 )
    return 1;

  dataFile >> compressionValue;
  dataFile >> totalNumberOfFile;
 
  for ( G4int i = 1;i <= totalNumberOfFile;i++ )
    {
      dataFile >> nameOfFileBuffer;
      listOfFile.push_back( nameOfFileBuffer );
    }
			
  dataFile.close();
  return 0;
}

G4int DicomConfiguration::ReadG4File( G4String g4File )
{
  densityValue.clear();
	
  g4File = g4File + ".g4";
  std::ifstream readingG4FileHeader( g4File.c_str() );

  if ( readingG4FileHeader.good() != 1 )
    return 1;
		
  readingG4FileHeader >> totalRows >> totalColumns;
  readingG4FileHeader >> xPixelSpacing >> yPixelSpacing; 
  // X is horizontal, Y is vertical
  readingG4FileHeader >> sliceTickness;
  readingG4FileHeader >> sliceLocation;
  readingG4FileHeader >> compressionUsed;

  G4double densityValueBuffer = 0;
  while ( readingG4FileHeader >> densityValueBuffer )
    {
      densityValue.push_back( densityValueBuffer );
    }
	
  readingG4FileHeader.close();
  return 0;
}

G4double DicomConfiguration::GetDensityValue(G4int i)
{ 
  G4double value = 0.;
  if (i <0.) G4cout << "out of range in GetDensityValue()!"<<G4endl; 
  if ( i>=0) 
    {
     unsigned int j = i;
     if(j >= densityValue.size() )
    {
      // Throw exception, return dummy, cerr error message...
      G4cout << "out of range in GetDensityValue()!"<<G4endl;
    }
  else
    {
      value = densityValue[i];
    }
    }
  return value;
}
