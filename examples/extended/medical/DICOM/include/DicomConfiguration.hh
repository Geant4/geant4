//   $tigre.2@sympatico.ca
//   08/05/03

//*******************************************************
//
// DicomConfiguration.hh :
//	- Handling of the header of *.g4 files
//	- Reading <Data.dat> file
//	- Used each we need to know information about the
//	  picture we are analysing
//	- Functions are in DicomConfiguration.cc
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



#ifndef DicomConfiguration_h
#define DicomConfiguration_h 1

// Standard library
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>

using namespace std;

class DicomConfiguration
{
public:

	// This function reads <Data.dat>, return 0/1 if successfull or not
	int Read_DataFile();
	
	// Public variables used in Read_DataFile(), those are the variables that contains informations
	int CompressionValue;
	int TotalNumberOfFile;
	vector<string> ListOfFile;
	
	// This function reads one *.g4 given by char* g4File, return 0/1 if successfull or not
	int Read_g4File(string g4File);
	
	// Public variables used in Read_g4File(), those are the variables that contains informations
	int TotalRows, TotalColumns;
	double  X_PixelSpacing, Y_PixelSpacing;
	double  SliceTickness;
	double  SliceLocation;
	int CompressionUsed;
	vector<double> DensityValue;
	
private:

	// Private variables used in Read_DataFile()
	string NameOfFileBuffer;
	
}
;//$   vihut@phy.ulaval.ca

#endif

