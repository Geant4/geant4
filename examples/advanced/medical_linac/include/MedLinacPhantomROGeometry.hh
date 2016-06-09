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
// $Id: MedLinacPhantomROGeometry.hh,v 1.3 2006/06/29 16:03:53 gunter Exp $
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
