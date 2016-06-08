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
// $Id: B08ImportanceDetectorConstruction.hh,v 1.1 2002/06/04 11:14:51 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef B08ImportanceDetectorConstruction_hh 
#define B08ImportanceDetectorConstruction_hh  B08ImportanceDetectorConstruction_hh 

#include "globals.hh"

class G4VPhysicalVolume;
class G4VIStore;

class B08ImportanceDetectorConstruction
{
public:
  B08ImportanceDetectorConstruction();
  ~B08ImportanceDetectorConstruction();

  G4VIStore* GetIStore();
  G4VPhysicalVolume *GetWorldVolume();

private:
  void Construct();
  G4VIStore *fnIStore;
  G4VPhysicalVolume *fWorldVolume;

  
};

#endif
