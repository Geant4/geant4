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
// $Id: B02DetectorConstruction.hh,v 1.3 2002-11-08 14:47:42 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef B02DetectorConstruction_hh
#define B02DetectorConstruction_hh B02DetectorConstruction_hh

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4IStore;

class B02DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B02DetectorConstruction();
  ~B02DetectorConstruction();
  
  G4VPhysicalVolume* Construct();

};

#endif
