
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

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "globals.hh"
#include "g4std/vector"
using namespace std;

class DicomConfiguration
{
public:
  DicomConfiguration(){;}
  ~DicomConfiguration(){;}
public:
  // This function reads <Data.dat>, return 0/1 if successfull or not
  G4int ReadDataFile();
  G4int ReadG4File(string g4File);
  G4int GetCompressionValue(){return CompressionValue;}
  G4int GetTotalNumberOfFile(){return TotalNumberOfFile;}
  G4std::vector<G4String> GetListOfFile(){return ListOfFile;}
  G4int GetTotalRows(){return TotalRows;}
  G4int GetTotalColumns(){return TotalColumns;}  
  G4double GetXPixelSpacing(){return  X_PixelSpacing;}
  G4double GetYPixelSpacing(){return  Y_PixelSpacing;}
  G4double GetSliceThickness(){return SliceTickness;} 
  G4double GetSliceLocation(){return  SliceLocation;} 
  G4int IsCompressionUsed(){return CompressionUsed;}
 
private:
	
  G4int CompressionValue;
  G4int TotalNumberOfFile;
  G4std::vector<G4String> ListOfFile;
  G4int TotalRows;
  G4int TotalColumns;
  G4double X_PixelSpacing;
  G4double Y_PixelSpacing;
  G4double SliceTickness;
  G4double SliceLocation; 
  G4int CompressionUsed;
  G4String NameOfFileBuffer;
	
public:
  G4std::vector<G4double> DensityValue;	
};

#endif

