//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------------
// $Id: MicrobeamPhantomConfiguration.hh,v 1.4 2006-06-28 13:43:05 gunter Exp $
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
   G4float GetCytoplasmMass() {return cytoplasmMass;}  
   G4float GetNucleusMass()   {return nucleusMass;}  

   G4ThreeVector GetVoxelThreeVector(G4int i) {return voxelThreeVector[i];}

private:
   
   static G4int phantomTotalPixels;
   static G4int nucleusTotalPixels;
   static G4int cytoplasmTotalPixels;
   static G4float dx;
   static G4float dy;
   static G4float dz;
   static G4float nucleusMass;
   static G4float cytoplasmMass;
   
   G4ThreeVector * voxelThreeVector;

};
#endif

