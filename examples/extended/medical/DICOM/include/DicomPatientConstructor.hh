#ifndef DicomPatientConstructor_h
#define DicomPatientConstructor_h 1

#include "globals.hh"

class DicomPatientConstructor
{
public:

  DicomPatientConstructor() {;}

  ~DicomPatientConstructor() {;}

  G4int FindingNbOfVoxels(G4double maxDensity, G4double minDensity);

  // Functions to use ROI (region of interest), contour usually drawn by the
  // physician to identify tumor volume and organ at risk
  void readContour();

  G4bool isWithin(G4double x, G4double y, G4double z);

private:  

  G4double contoursX[100][100];
  G4double contoursY[100][100];
  G4double contoursZ[100][100];
  G4int maxCurve;

  //char  maxBuf[300];
  //char compressionBuf[300];
  //char name[300];
  //G4int compression;
  //G4double pixelSpacingX,pixelSpacingY;
  //G4double sliceTickness;
  //G4double sliceLocation;
  //char rowsBuf[300],columnsBuf[300];
  //char pixelSpacingXBuf[300],pixelSpacingYBuf[300];
  //char sliceTicknessBuf[300];
  //char sliceLocationBuf[300];
  //char fullname[300];
  
  //FILE* readConf;
  //G4int flag_contours;
  //G4int lenc,lenr;
  //char densityBuf[300];
  //G4std::vector<G4double> density;

};
#endif
