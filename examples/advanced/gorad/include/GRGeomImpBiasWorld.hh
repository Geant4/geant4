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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRGeomImpBiasWorld.hh
//   Header file of a parallel world class that defines the geometry
//   of the geometry improtance biasing.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRGeomImpBiasWorld_H
#define GRGeomImpBiasWorld_H 1

#include "G4VUserParallelWorld.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class GRDetectorConstruction;
class G4VPhysicalVolume;
class G4Region;
class GRGeomImpBiasWorldStateNotifier;

class GRGeomImpBiasWorld : public G4VUserParallelWorld
{
  friend class GRGeomImpBiasWorldStateNotifier;

  public:
    GRGeomImpBiasWorld(G4String&,GRDetectorConstruction*);
    virtual ~GRGeomImpBiasWorld();
    virtual void Construct();
    virtual void ConstructSD();

  private:
    GRDetectorConstruction* detector;
    G4VPhysicalVolume* fWorld = nullptr;
    G4Region* biasingRegion = nullptr;
    GRGeomImpBiasWorldStateNotifier* stateNotifier;
    G4bool constructed = false;

  private:
    void Update();
};

#include "G4VStateDependent.hh"
#include "G4ApplicationState.hh"

class GRGeomImpBiasWorldStateNotifier : public G4VStateDependent
{
  public:
    GRGeomImpBiasWorldStateNotifier(GRGeomImpBiasWorld* bw)
    : biasWorld(bw)
    {;}
    virtual ~GRGeomImpBiasWorldStateNotifier()
    {;}
    virtual G4bool Notify(G4ApplicationState newState)
    {
      if(previousState == G4State_Idle && newState == G4State_GeomClosed)
      { biasWorld->Update(); }
      previousState = newState;
      return true;
    }
  private:
    GRGeomImpBiasWorld* biasWorld;
    G4ApplicationState previousState = G4State_PreInit;
};

#endif
