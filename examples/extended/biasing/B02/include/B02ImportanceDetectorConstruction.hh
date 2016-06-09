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
//
// $Id: B02ImportanceDetectorConstruction.hh,v 1.7 2007-06-22 13:38:55 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef B02ImportanceDetectorConstruction_hh 
#define B02ImportanceDetectorConstruction_hh  B02ImportanceDetectorConstruction_hh 

#include "globals.hh"
#include <map>
#include <vector>
#include "G4GeometryCell.hh"
#include "B02PVolumeStore.hh"

#include "G4VUserParallelWorld.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;



class B02ImportanceDetectorConstruction : public G4VUserParallelWorld
{
public:
  B02ImportanceDetectorConstruction(G4String worldName);
  ~B02ImportanceDetectorConstruction();

  const G4VPhysicalVolume &GetPhysicalVolumeByName(const G4String& name) const;
  G4VPhysicalVolume &GetWorldVolumeAddress() const;
  G4String ListPhysNamesAsG4String();
  G4String GetCellName(G4int i);
  G4GeometryCell GetGeometryCell(G4int i);

  G4VPhysicalVolume* GetWorldVolume();

  void SetSensitive();

private:
  void Construct();
  B02PVolumeStore fPVolumeStore;

  //  std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;
  std::vector< G4LogicalVolume * > fLogicalVolumeVector;

  //  G4VPhysicalVolume *fWorldVolume;

  G4VPhysicalVolume* ghostWorld;

};

#endif
