//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------------
// $Id: MicrobeamPhantomConfiguration.hh,v 1.5 2006-06-29 16:05:07 gunter Exp $
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

