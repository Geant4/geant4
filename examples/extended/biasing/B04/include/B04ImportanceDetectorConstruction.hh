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
// $Id: B04ImportanceDetectorConstruction.hh,v 1.2 2002/04/19 10:54:29 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef B04ImportanceDetectorConstruction_hh 
#define B04ImportanceDetectorConstruction_hh  B04ImportanceDetectorConstruction_hh 

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4VIStore;

class B04ImportanceDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B04ImportanceDetectorConstruction();
  ~B04ImportanceDetectorConstruction();

  G4VPhysicalVolume* Construct();
  G4VIStore* GetIStore();

private:
  G4VIStore *fIStore;
  
};

#endif
