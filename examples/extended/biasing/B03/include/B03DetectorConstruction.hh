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
/// \file biasing/B03/include/B03DetectorConstruction.hh
/// \brief Definition of the B03DetectorConstruction class
//
//
//

#ifndef B03DetectorConstruction_hh
#define B03DetectorConstruction_hh B03DetectorConstruction_hh

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4IStore;

class B03DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B03DetectorConstruction();
  ~B03DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

  G4VPhysicalVolume* GetWorldVolume();
  G4VPhysicalVolume& GetWorldVolumeAddress() const;

  //  G4String GetCellName(G4int i);

  void SetSensitive();

  //  virtual void ConstructSDandField();

private:

  G4VPhysicalVolume* fWorldVolume;

};

#endif
