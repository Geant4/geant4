#ifndef DicomPatientConstructor_h
#define DicomPatientConstructor_h 1

#include "globals.hh"
#include <stdio.h>
#include <vector>
using namespace std;

class DicomPatientConstructor
{
public:
  DicomPatientConstructor(){;}
  ~DicomPatientConstructor(){;}
  G4int FindingNbOfVoxels(G4double MaxDensity , G4double MinDensity);

 // Functions to use ROI (region of interest), contour usually drawn by the
  // physician to identify tumor volume and organ at risk
  void readContour();
  G4bool isWithin(G4double,G4double,G4double);
  
  double ContoursX[100][100];
  double ContoursY[100][100];
  double ContoursZ[100][100];
  int MaxCurve;  
  int max;
  char maxbuf[300];
  int compression;
  char compressionbuf[300];
  char name[300];
  int columns,rows;
  double pixel_spacing_X,pixel_spacing_Y;
  double SliceTickness;
  double SliceLocation;
  char rowsbuf[300],columnsbuf[300];
  char pixel_spacing_Xbuf[300],pixel_spacing_Ybuf[300];
  char SliceTicknessbuf[300];
  char SliceLocationbuf[300];
  char fullname[300];
  FILE* readData;
  FILE* readConf;
 int flag_contours;
  int lenc,lenr;
  char Densitybuf[300];
  vector<double> Density;

};
#endif
