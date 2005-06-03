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
// $Id: HadrontherapyPhantomROGeometry.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

//The phantom is devided in voxels. the dimension of the voxel is 1 mm
//
//
#ifndef HadrontherapyPhantomROGeometry_h
#define HadrontherapyPhantomROGeometry_h 

#include "G4VReadOutGeometry.hh"

class HadrontherapyPhantomROGeometry : public G4VReadOutGeometry
{
public:
  HadrontherapyPhantomROGeometry(G4String aString,
				 G4double phantomDimX,
				 G4double phantomDimY,
				 G4double phantomDimZ,
				 G4int numberOfVoxelsX,
				 G4int numberOfVoxelsY,
				 G4int numberOfVoxelsZ);

  ~HadrontherapyPhantomROGeometry();

private:
  G4VPhysicalVolume* Build();

private:  
  const G4double phantomSizeX;
  const G4double phantomSizeY; 
  const G4double phantomSizeZ;

  const G4int numberOfVoxelsAlongX;
  const G4int numberOfVoxelsAlongY; 
  const G4int numberOfVoxelsAlongZ; 
  
  G4VPhysicalVolume *ROPhantomZDivisionPhys;
};
#endif
