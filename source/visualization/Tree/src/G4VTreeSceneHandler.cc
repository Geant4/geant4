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
// $Id: G4VTreeSceneHandler.cc,v 1.7 2002-12-11 16:08:55 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4VTreeGraphicsScene (was in modeling).

#include "G4VTreeSceneHandler.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ModelingParameters.hh"

G4int G4VTreeSceneHandler::fSceneIdCount = 0;
// Counter for Tree scene handlers.

G4int G4VTreeSceneHandler::fSceneCount = 0;
// No. of extanct scene handlers.

G4VTreeSceneHandler::G4VTreeSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name),
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0),
  fpCurrentObjectTransformation (0)
{
  fSceneCount++;
}

G4VTreeSceneHandler::~G4VTreeSceneHandler () {}

void G4VTreeSceneHandler::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4VTreeSceneHandler::BeginModeling() {
  G4VSceneHandler::BeginModeling();  // Required: see G4VSceneHandler.hh.
}

void G4VTreeSceneHandler::EndModeling() {
  G4VSceneHandler::EndModeling();  // Required: see G4VSceneHandler.hh.
}
