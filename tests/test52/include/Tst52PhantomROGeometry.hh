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
//
// $Id: Tst52PhantomROGeometry.hh,v 1.1.2.1 2007-12-10 16:33:48 gunter Exp $
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
