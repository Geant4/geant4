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

// Definitions
#include "DicomConfiguration.hh"


using namespace::std;

G4int DicomConfiguration::ReadDataFile()
{
	ifstream DataFile( "Data.dat" );

	if ( DataFile.good() != 1 )
		return 1;

	DataFile >> CompressionValue;
	DataFile >> TotalNumberOfFile;

	for(int i=1;i<=TotalNumberOfFile;i++)
	{
		DataFile >> NameOfFileBuffer;
		ListOfFile.push_back( NameOfFileBuffer );
	}
			
	DataFile.close();
	return 0;
}

G4int DicomConfiguration::ReadG4File( string g4File )
{

	DensityValue.clear();
	
	g4File=g4File + ".g4";
	ifstream Reading_g4FileHeader( g4File.c_str() );

	if ( Reading_g4FileHeader.good() != 1 )
		return 1;
		
	Reading_g4FileHeader >> TotalRows >> TotalColumns;
	Reading_g4FileHeader >> X_PixelSpacing >> Y_PixelSpacing; // X is horizontal, Y is vertical
	Reading_g4FileHeader >> SliceTickness;
	Reading_g4FileHeader >> SliceLocation;
	Reading_g4FileHeader >> CompressionUsed;

	double DensityValueBuffer=0;
	while ( Reading_g4FileHeader >> DensityValueBuffer )
	{
		DensityValue.push_back( DensityValueBuffer );
	}
	
	Reading_g4FileHeader.close();
	return 0;
}
