
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

#include "globals.hh"
#include "g4std/vector"

class DicomConfiguration
{
public:

  DicomConfiguration() { }
  ~DicomConfiguration() { }

  // This function reads <Data.dat>, return 0/1 if successfull or not
  G4bool ReadDataFile();

  G4int ReadG4File(G4String g4File);

  G4int GetCompressionValue() {return compressionValue;}

  G4int GetTotalNumberOfFile() {return totalNumberOfFile;}

  G4std::vector<G4String> GetListOfFile() {return listOfFile;}

  G4int GetTotalRows() {return totalRows;}

  G4int GetTotalColumns() {return totalColumns;}  

  G4double GetXPixelSpacing() {return  xPixelSpacing;}

  G4double GetYPixelSpacing() {return  yPixelSpacing;}

  G4double GetSliceThickness() {return sliceTickness;} 

  G4double GetSliceLocation() {return  sliceLocation;} 

  G4int IsCompressionUsed() {return compressionUsed;}

  G4double GetDensityValue(G4int i);
 
private:
	
  G4int compressionValue;
  G4int totalNumberOfFile;
  G4std::vector<G4String> listOfFile;
  G4int totalRows;
  G4int totalColumns;
  G4double xPixelSpacing;
  G4double yPixelSpacing;
  G4double sliceTickness;
  G4double sliceLocation; 
  G4int compressionUsed;
  G4String nameOfFileBuffer;	
  G4std::vector<G4double> densityValue;	

};

#endif

