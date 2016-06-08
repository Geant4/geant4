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
// $Id: B07DetectorConstruction.hh,v 1.2 2002/04/19 10:54:34 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef B07DetectorConstruction_hh
#define B07DetectorConstruction_hh B07DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class B07DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B07DetectorConstruction();
  ~B07DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
};

#endif
