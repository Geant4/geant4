//   $tigre.2@sympatico.ca, louis.archambault@phy.ulaval.ca
//   03/10/02

//*******************************************************
//
// DicomHandler.hh :
//	- Handling of DICM images
//	- Transforming *.dcm to *.g4 ( pixels->density )
//	- Reading headers and pixels
//	- Transforming pixel to density and creating *.g4
//	  files
//	- Functions are in DicomHandler.cc
//
// The code was written by :
//	Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
// Base on previous code by :
//	Dragan Tubic <tdragan@gel.ulaval.ca>
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

#ifndef DicomHandler_h
#define DicomHandler_h 1

#include "globals.hh"

class DicomHandler
{
public:

  DicomHandler();

  ~DicomHandler() { }

  G4int readHeader(FILE *,char[300]);

  G4int readData(FILE *,char[300]); // note: always use readHeader 
                                    // before readData

  // use ImageMagick to display the image
  //G4int displayImage(char[500]);

  void checkFileFormat();

private:

  G4int compression, max;
  G4double pixel2density(G4int pixel);
  G4int rows;
  G4int columns;
  G4int bitAllocated;

  char  pixelSpacingX[300],pixelSpacingY[300];
  char sliceThickness[300];
  G4double sliceLocation;
};
#endif

