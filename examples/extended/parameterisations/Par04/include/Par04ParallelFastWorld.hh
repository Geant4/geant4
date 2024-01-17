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
#ifndef PAR04PARALLELFASTWORLD_HH
#define PAR04PARALLELFASTWORLD_HH

#include "Par04DetectorConstruction.hh"
#include "G4VUserParallelWorld.hh"
#include "globals.hh"

class Par04ParallelMessenger;
class Par04ParallelFullWorld;
class G4LogicalVolume;
class G4VPhysicalVolume;

class Par04ParallelFastWorld : public G4VUserParallelWorld
{
public:
  Par04ParallelFastWorld(G4String aWorldName,
                         const Par04DetectorConstruction* aMassDetector,
                         const Par04ParallelFullWorld* aParallelFull);
  ~Par04ParallelFastWorld();
  
  virtual void Construct() final;
  virtual void ConstructSD() final;
  
  void Print();

private:
  ///  Messenger that allows to modify geometry
  const Par04DetectorConstruction* fMassDetector;
  const Par04ParallelFullWorld* fParallelFull;
  std::vector<G4LogicalVolume*> fLogicalCell;
  G4int fNbOfLayers = 1;
  G4int fNbOfSlices = 1;
  G4int fNbOfRows = 1;
  G4double fLayerThickness = 0;
};

#endif
