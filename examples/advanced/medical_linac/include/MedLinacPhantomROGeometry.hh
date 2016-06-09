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
// $Id: MedLinacPhantomROGeometry.hh,v 1.2 2004/04/02 17:48:41 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
//The phantom is devided in voxels. the dimension of the voxel is 1mm
//
#ifndef MedLinacPhantomROGeometry_h
#define MedLinacPhantomROGeometry_h 

#include "G4VReadOutGeometry.hh"

class MedLinacPhantomROGeometry : public G4VReadOutGeometry
{
public:
  MedLinacPhantomROGeometry(G4String aString,
			  G4double phantomDim_x,
			  G4double phantomDim_y,
			  G4double phantomDim_z,
			  G4int numberOfVoxelsX,
			  G4int numberOfVoxelsY,
			  G4int numberOfVoxelsZ);
  ~MedLinacPhantomROGeometry();
private:
  G4VPhysicalVolume* Build();

private:
  const G4double PhantomDimensionX;
  const G4double PhantomDimensionY;
  const G4double PhantomDimensionZ;
  const G4int NumberOfVoxelsAlongX;
  const G4int NumberOfVoxelsAlongY;
  const G4int NumberOfVoxelsAlongZ; 
  G4VPhysicalVolume *ROPhantomYDivisionPhys;
};
#endif
