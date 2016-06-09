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
// $Id: RemSimROGeometry.hh,v 1.5 2004/05/22 12:57:05 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//    ************************************
//    *                                  *
//    *            RemSimROGeometry.hh   *
//    *                                  *
//    ************************************
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//

#ifndef RemSimROGeometry_h
#define RemSimROGeometry_h 

#include "G4VReadOutGeometry.hh"

class RemSimROGeometry : public G4VReadOutGeometry
{
public:
  RemSimROGeometry( G4double phantomDimX,
                    G4double phantomDimY,
		    G4double phantomDimZ,
		    G4int numberOfVoxelsZ,
                    G4double trans);

  ~RemSimROGeometry();

private:
  G4VPhysicalVolume* Build();

private:
  const G4double astronautDimensionX;
  const G4double astronautDimensionY;
  const G4double astronautDimensionZ;
  const G4int numberOfVoxelsAlongZ; 
  const G4double translation; 
  G4VPhysicalVolume *ROAstronautZDivisionPhys; 
};
#endif
