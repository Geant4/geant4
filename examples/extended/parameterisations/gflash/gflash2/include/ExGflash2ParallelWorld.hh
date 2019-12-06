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
//
#ifndef ExGflash2ParallelWorld_h
#define ExGflash2ParallelWorld_h 1

#include "G4VUserParallelWorld.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Region;

class GFlashHomoShowerParameterisation;
class GFlashShowerModel;
class GFlashHitMaker;
class GFlashParticleBounds;

class ExGflash2ParallelWorld : public G4VUserParallelWorld
{
public:
  ExGflash2ParallelWorld(G4String aWorldName);
  ~ExGflash2ParallelWorld();
  
  virtual void Construct() final;
  virtual void ConstructSD() final;
  
  
private:
  G4Region*           fRegion;

  static G4ThreadLocal GFlashShowerModel* fFastShowerModel; 
  static G4ThreadLocal GFlashHomoShowerParameterisation* fParameterisation;
  static G4ThreadLocal GFlashParticleBounds* fParticleBounds;
  static G4ThreadLocal GFlashHitMaker* fHitMaker;
};

#endif




















