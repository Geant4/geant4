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
//
// $Id: Tst52PhantomROGeometry.hh,v 1.1 2007-04-12 12:00:17 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    ************************************
//    *                                  *
//    *    Tst52PhantomROGeometry.hh   *
//    *                                  *
//    ************************************
//Author: Susanna Guatelli (guatelli@ge.infn.it)
//The phantom is devided in voxels. the dimension of the voxel is 1mm
//
#ifndef Tst52PhantomROGeometry_h
#define Tst52PhantomROGeometry_h 

#include "G4VReadOutGeometry.hh"

class Tst52PhantomROGeometry : public G4VReadOutGeometry
{
public:
  Tst52PhantomROGeometry(G4String aString);

  ~Tst52PhantomROGeometry();

private:
  G4VPhysicalVolume* Build();

public:
  void SetROParameter(G4double phantomDimX,
		      G4double phantomDimY,
		      G4double phantomDimZ,
		      G4int numberOfVoxelsZ);
private:
  G4double phantomDimensionX;
  G4double phantomDimensionY;
  G4double phantomDimensionZ;

  G4int numberOfVoxelsAlongZ; 
  G4VPhysicalVolume *ROPhantomZDivisionPhys;
};
#endif
