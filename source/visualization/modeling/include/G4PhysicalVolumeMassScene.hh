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
// $Id: G4PhysicalVolumeMassScene.hh 99076 2016-09-01 12:40:47Z gcosmo $
//
// 
// John Allison  12th September 2004

// Class Description:
//
// Calculates the mass of a geometry tree taking into account daughters
// up to the depth specified in the G4PhysicalVolumeModel.  Culling is
// ignored so that all volumes are seen.
//
// Do not use this for a "parallel world" for which materials are not
// defined.  Use only for the material world.
//
// The calculation is quite tricky, since it involves subtracting the
// mass of that part of the mother that is occupied by each daughter and
// then adding the mass of the daughter, and so on down the heirarchy.
//
// Usage for a given G4PhysicalVolumeModel* pvModel:
//   G4PhysicalVolumeMassScene massScene(pvModel);
//   pvModel->DescribeYourselfTo (massScene);
//   G4double volume = massScene.GetVolume();
//   G4double mass = massScene.GetMass();
//   massScene.Reset();
// See, for example, G4ASCIITreeSceneHandler::EndModeling().

#ifndef G4PHYSICALVOLUMEMASSSCENE_HH
#define G4PHYSICALVOLUMEMASSSCENE_HH

#include "G4PseudoScene.hh"

#include <deque>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PhysicalVolumeModel;
class G4Material;

class G4PhysicalVolumeMassScene: public G4PseudoScene
{

public:
  G4PhysicalVolumeMassScene (G4PhysicalVolumeModel*);
  virtual ~G4PhysicalVolumeMassScene ();

public: // With description

  G4double GetVolume () const {return fVolume;}
  // Overall volume.

  G4double GetMass () const {return fMass;}
  // Mass of whole tree, i.e., accounting for all daughters.

  void Reset ();
  // Reset for subsequent re-use.

private:
  void ProcessVolume (const G4VSolid&);
  G4PhysicalVolumeModel* fpPVModel;
  G4double fVolume;
  G4double fMass;
  G4VPhysicalVolume* fpLastPV;
  G4int fPVPCount;
  G4int fLastDepth;
  G4double fLastDensity;
  std::deque<G4double> fDensityStack;
};

#endif
