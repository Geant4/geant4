#ifndef DicomPatientConstructor_h
#define DicomPatientConstructor_h 1

#include "globals.hh"
#include "g4std/vector"

class DicomPatientConstructor
{
public:

  DicomPatientConstructor() {;}
  ~DicomPatientConstructor() {;}

  G4int FindingNbOfVoxels(G4double maxDensity, G4double minDensity);

  // Functions to use ROI (region of interest), contour usually drawn by the
  // physician to identify tumor volume and organ at risk 
  // under development ...
  // void readContour();
  //G4bool isWithin(G4double,G4double,G4double);

private:  
  G4double pixelSpacingX;
  G4double pixelSpacingY;
  G4double sliceThickness;
  G4double sliceLocation;
  G4int lenc,lenr;
};
#endif
