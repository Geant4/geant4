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
/// \file biasing/B01/include/B01DetectorConstruction.hh
/// \brief Definition of the B01DetectorConstruction class
//
//
//

#ifndef B01DetectorConstruction_hh
#define B01DetectorConstruction_hh B01DetectorConstruction_hh

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <vector>
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VIStore;
class G4VWeightWindowStore;

class B01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B01DetectorConstruction();
  ~B01DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

  G4VIStore* CreateImportanceStore();
    // create an importance store, caller is responsible for deleting it

  G4VWeightWindowStore *CreateWeightWindowStore();
    // create an weight window  store, caller is responsible for 
    // deleting it

  G4String GetCellName(G4int i);

  G4VPhysicalVolume* GetWorldVolume();

  void SetSensitive();

  virtual void ConstructSDandField();

private:
  std::vector< G4LogicalVolume * > fLogicalVolumeVector;
  std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;

  G4VPhysicalVolume* fWorldVolume;

};

#endif
