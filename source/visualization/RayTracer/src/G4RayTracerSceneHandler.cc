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
// $Id: G4RayTracerSceneHandler.cc,v 1.3 2005/06/02 17:43:46 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $

#include "G4RayTracerSceneHandler.hh"

G4RayTracerSceneHandler::G4RayTracerSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name)
{}

G4RayTracerSceneHandler::~G4RayTracerSceneHandler()
{}

G4int G4RayTracerSceneHandler::fSceneIdCount = 0;
