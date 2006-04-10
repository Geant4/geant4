// -------------------------------------------------------------------
// $Id: MicrobeamPhantomConfiguration.hh,v 1.2 2006-04-10 14:47:31 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamPhantomConfiguration_h
#define MicrobeamPhantomDicomConfiguration_h 1

#include "G4ThreeVector.hh"

class MicrobeamPhantomConfiguration
{
public:

   MicrobeamPhantomConfiguration();
  ~MicrobeamPhantomConfiguration();

   G4int   Initialize();
   G4int   GetPhantomTotalPixels()   {return phantomTotalPixels;}  
   G4int   GetNucleusTotalPixels()   {return nucleusTotalPixels;}  
   G4int   GetCytoplasmTotalPixels() {return cytoplasmTotalPixels;}  
   G4float GetPixelSizeX() {return dx;}  
   G4float GetPixelSizeY() {return dy;}  
   G4float GetPixelSizeZ() {return dz;}  

   G4ThreeVector GetVoxelThreeVector(G4int i) {return voxelThreeVector[i];}

private:
   
   static G4int phantomTotalPixels;
   static G4int nucleusTotalPixels;
   static G4int cytoplasmTotalPixels;
   static G4float dx;
   static G4float dy;
   static G4float dz;
   
   G4ThreeVector * voxelThreeVector;

};
#endif

