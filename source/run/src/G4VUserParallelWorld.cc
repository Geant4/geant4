//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VUserParallelWorld.cc,v 1.1 2006-04-26 15:24:24 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VUserParallelWorld.hh"

#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

G4VUserParallelWorld::G4VUserParallelWorld(G4String worldName)
{ fWorldName = worldName; }

G4VUserParallelWorld::~G4VUserParallelWorld()
{ ; }

G4VPhysicalVolume* G4VUserParallelWorld::GetWorld()
{
  G4VPhysicalVolume* pWorld
       = G4TransportationManager::GetTransportationManager()
         ->GetParallelWorld(fWorldName);
  pWorld->SetName(fWorldName);
  return pWorld;
}


