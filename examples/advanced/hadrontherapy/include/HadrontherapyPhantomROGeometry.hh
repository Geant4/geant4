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
//The phantom is devided in voxels. the dimension of the voxel is 1mm
//
#ifndef HadrontherapyPhantomROGeometry_h
#define HadrontherapyPhantomROGeometry_h 

#include "G4VReadOutGeometry.hh"

class HadrontherapyPhantomROGeometry : public G4VReadOutGeometry
{
public:
  HadrontherapyPhantomROGeometry(G4String aString,
			  G4double phantomDimX,
			  G4double phantomDimZ,
			  G4int numberOfVoxelsX,
			  G4int numberOfVoxelsZ);
  ~HadrontherapyPhantomROGeometry();
private:
  G4VPhysicalVolume* Build();

private:
  const G4double phantomDimensionX;
  const G4double phantomDimensionZ;


  const G4int numberOfVoxelsAlongX;
  const G4int numberOfVoxelsAlongZ; 
  G4VPhysicalVolume *ROPhantomYDivisionPhys;
};
#endif
