
//
// DicomConfiguration.cc :
//	- Handling of the header of *.g4 files
// 	- Reading <Data.dat> file
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
#include "DicomConfiguration.hh"
#include "g4std/fstream"
#include "g4std/vector"
G4bool DicomConfiguration::ReadDataFile()
{
  G4std::ifstream dataFile( "Data.dat" );

  if ( dataFile.good() != 1 )
    return 1;

  dataFile >> compressionValue;
  dataFile >> totalNumberOfFile;

  for (G4int i=1;i<=totalNumberOfFile;i++)
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
  G4std::ifstream readingG4FileHeader( g4File.c_str() );

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
  if (i < 0 || i >= densityValue.size() )
    {
      // Throw exception, return dummy, cerr error message...
      G4cout << "out of range in GetDensityValue()!"<<G4endl;
    }
  else
    {
      value = densityValue[i];
    }
 
  return value;
}
