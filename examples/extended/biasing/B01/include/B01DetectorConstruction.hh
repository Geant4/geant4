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
// $Id: B01DetectorConstruction.hh,v 1.3 2002-09-17 13:59:15 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef B01DetectorConstruction_hh
#define B01DetectorConstruction_hh B01DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class B01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B01DetectorConstruction(G4VPhysicalVolume &worldvol);
  ~B01DetectorConstruction();
  
  G4VPhysicalVolume* Construct();

private:
  G4VPhysicalVolume* fWorldVolume;
};

#endif
