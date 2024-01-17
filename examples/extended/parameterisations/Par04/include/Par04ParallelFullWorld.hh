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
#ifndef PAR04PARALLELFULLWORLD_HH
#define PAR04PARALLELFULLWORLD_HH

#include "Par04DetectorConstruction.hh"
#include "G4VUserParallelWorld.hh"
#include "globals.hh"

class Par04ParallelMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

class Par04ParallelFullWorld : public G4VUserParallelWorld
{
public:
  Par04ParallelFullWorld(G4String aWorldName, const Par04DetectorConstruction* aMassDetector);
  ~Par04ParallelFullWorld();
  
  virtual void Construct() final;
  virtual void ConstructSD() final;

  
  ///  Set number of slices
  inline void SetNbOfSlices(G4int aNumber) { fNbOfSlices = aNumber; };
  ///  Get number of slices
  inline G4int GetNbOfSlices() const { return fNbOfSlices; };
  ///  Set number of rows
  inline void SetNbOfRows(G4int aNumber) { fNbOfRows = aNumber; };
  ///  Get number of rows
  inline G4int GetNbOfRows() const { return fNbOfRows; };
  ///  Get number of layers
  inline G4int GetNbOfLayers() const { return fNbOfLayers; };
  
  void Print();

private:
  ///  Messenger that allows to modify geometry
  Par04ParallelMessenger* fParallelMessenger = nullptr;
  const Par04DetectorConstruction* fMassDetector;
  std::vector<G4LogicalVolume*> fLogicalCell;
  G4int fNbOfLayers = 1;
  G4int fNbOfSlices = 1;
  G4int fNbOfRows = 1;
  G4double fLayerThickness = 0;
};

#endif
