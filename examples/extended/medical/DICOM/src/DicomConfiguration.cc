//   $tigre.2@sympatico.ca
//   08/05/03

//*******************************************************
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


#include "DicomConfiguration.hh"
#include "g4std/fstream"

G4int DicomConfiguration::ReadDataFile()
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

G4int DicomConfiguration::ReadG4File( const G4String& g4File )
{
  densityValue.clear();
	
  G4String file = g4File + ".g4";
  G4std::ifstream readingG4FileHeader( file.c_str() );

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
  // ---- MGP ---- Decide logic for dummy or exception
  G4double value = 0.;
  if (i < 0 || i >= densityValue.size() )
    {
      // Throw exception, return dummy, cerr error message...
    }
  else
    {
      value = densityValue[i];
    }
 
  return value;
}
