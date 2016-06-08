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
// $Id: B06ImportanceDetectorConstruction.hh,v 1.3 2002/04/19 10:54:31 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef B06ImportanceDetectorConstruction_hh 
#define B06ImportanceDetectorConstruction_hh  B06ImportanceDetectorConstruction_hh 

#include "globals.hh"

class G4VPhysicalVolume;
class G4VIStore;

class B06ImportanceDetectorConstruction
{
public:
  B06ImportanceDetectorConstruction();
  ~B06ImportanceDetectorConstruction();

  G4VIStore* GetIStore();

private:
  void Construct();
  G4VIStore *fIStore;
  
};

#endif
