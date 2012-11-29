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
/// \file biasing/B02/include/B02PSScoringDetectorConstruction.hh
/// \brief Definition of the B02PSScoringDetectorConstruction class
//
//
// $Id$
//

#ifndef B02PSScoringDetectorConstruction_hh 
#define B02PSScoringDetectorConstruction_hh  B02PSScoringDetectorConstruction_hh 

#include "globals.hh"
#include <map>
//#include "G4GeometryCell.hh"
//#include "B02PVolumeStore.hh"

#include "G4VUserParallelWorld.hh"

#include <vector> 
class G4VPhysicalVolume;
class G4LogicalVolume;


class B02PSScoringDetectorConstruction: public G4VUserParallelWorld
{
public:
  B02PSScoringDetectorConstruction(G4String worldName);
  ~B02PSScoringDetectorConstruction();

  //const G4VPhysicalVolume &GetPhysicalVolumeByName(const G4String& name) const;
  //  G4VPhysicalVolume *GetWorldVolume() const;
  //G4String ListPhysNamesAsG4String();
  G4String GetCellName(G4int i);
  //G4GeometryCell GetGeometryCell(G4int i);

  G4VPhysicalVolume* GetWorldVolume();
  G4VPhysicalVolume& GetWorldVolumeAddress() const;

  void SetSensitive();

private:
  void Construct();
  //B02PVolumeStore fPVolumeStore;
  G4VPhysicalVolume *fWorldVolume;

  G4VPhysicalVolume* ghostWorld;

  std::vector< G4LogicalVolume* > fLogicalVolumeVector;

};

#endif
