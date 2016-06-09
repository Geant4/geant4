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
// $Id: B01DetectorConstruction.hh,v 1.6 2003/08/25 12:32:29 dressel Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef B01DetectorConstruction_hh
#define B01DetectorConstruction_hh B01DetectorConstruction_hh

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <vector>
class G4VPhysicalVolume;
class G4VIStore;
class G4VWeightWindowStore;

class B01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B01DetectorConstruction();
  ~B01DetectorConstruction();
  
  G4VPhysicalVolume* Construct();

  G4VIStore* CreateImportanceStore();
    // create an importance store, caller is responsible for deleting it

  G4VWeightWindowStore *CreateWeightWindowStore();
    // create an weight window  store, caller is responsible for 
    // deleting it

  G4String GetCellName(G4int i);
private:
  std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;

};

#endif
