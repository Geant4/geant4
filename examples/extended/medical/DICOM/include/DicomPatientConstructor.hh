#ifndef DicomPatientConstructor_h
#define DicomPatientConstructor_h 1

#include "globals.hh"
#include "g4std/vector"
#include <stdio.h>

class DicomPatientConstructor
{
public:
  DicomPatientConstructor(){;}
  ~DicomPatientConstructor(){;}

  G4int FindingNbOfVoxels(G4double MaxDensity, G4double MinDensity);

 // Functions to use ROI (region of interest), contour usually drawn by the
  // physician to identify tumor volume and organ at risk
  void readContour();

  G4bool isWithin(G4double,G4double,G4double);

private:  
  G4double ContoursX[100][100];
  G4double ContoursY[100][100];
  G4double ContoursZ[100][100];
  G4int MaxCurve;
  char  maxbuf[300];
  char compressionbuf[300];
  char name[300];
  G4int compression;
  G4double pixel_spacing_X,pixel_spacing_Y;
  G4double SliceTickness;
  G4double SliceLocation;
  char rowsbuf[300],columnsbuf[300];
  char pixel_spacing_Xbuf[300],pixel_spacing_Ybuf[300];
  char SliceTicknessbuf[300];
  char SliceLocationbuf[300];
  char fullname[300];
  FILE* readData;
  FILE* readConf;
  G4int flag_contours;
  G4int lenc,lenr;
  char Densitybuf[300];
  G4std::vector<G4double> Density;

};
#endif
