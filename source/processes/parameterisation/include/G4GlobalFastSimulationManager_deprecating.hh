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
// $Id: G4GlobalFastSimulationManager_deprecating.hh,v 1.1 2006-11-03 17:33:30 mverderi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
//  G4GlobalFastSimulationManager_deprecating.hh
//
//  Description:
//    Gathers code to be dropped @ next major release
//
//  History:
//    November 06: Marc Verderi
//
//---------------------------------------------------------------

#ifndef  G4GlobalFastSimulationManager_deprecating_hh
#define  G4GlobalFastSimulationManager_deprecating_hh

#include "globals.hh"
#include "G4FastSimulationVector.hh"
#include "G4FlavoredParallelWorld.hh"
#include "G4StateManager.hh"

class G4GlobalFastSimulationManager;

class G4GlobalFastSimulationManager_deprecating {
public:

  ~G4GlobalFastSimulationManager_deprecating(); 

  void FastSimulationNeedsToBeClosed();

public:
  void CloseFastSimulation();
  G4VFlavoredParallelWorld* GetFlavoredWorldForThis(G4ParticleDefinition *);
  G4bool Notify(G4ApplicationState requestedState);

  G4GlobalFastSimulationManager_deprecating(G4GlobalFastSimulationManager* theGlobal);

private:
  G4bool fClosed;
  G4FastSimulationVector <G4FlavoredParallelWorld> NeededFlavoredWorlds;
  G4VPhysicalVolume* GiveMeAWorldVolumeClone();

  G4GlobalFastSimulationManager* _global;
};

#endif 
