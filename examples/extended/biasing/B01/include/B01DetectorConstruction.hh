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
// $Id: B01DetectorConstruction.hh,v 1.4 2002-10-22 14:09:02 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef B01DetectorConstruction_hh
#define B01DetectorConstruction_hh B01DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class B01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  explicit B01DetectorConstruction(G4VPhysicalVolume &worldvol);
  virtual ~B01DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

private:
  B01DetectorConstruction(const B01DetectorConstruction &);
  B01DetectorConstruction &operator=(const B01DetectorConstruction &);
  G4VPhysicalVolume* fWorldVolume;
};

#endif
