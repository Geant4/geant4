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
/// \file biasing/B02/include/B02DetectorConstruction.hh
/// \brief Definition of the B02DetectorConstruction class
//
//
// $Id$
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

  G4VPhysicalVolume* GetWorldVolume();
  G4VPhysicalVolume& GetWorldVolumeAddress() const;

  //  G4String GetCellName(G4int i);

  void SetSensitive();

private:

  G4VPhysicalVolume* pWorldVolume;


};

#endif
