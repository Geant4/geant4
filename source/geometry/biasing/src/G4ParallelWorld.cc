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
// $Id: G4ParallelWorld.cc,v 1.8 2006/06/29 18:17:33 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelWorld.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelWorld.hh"
#include "G4ParallelStepper.hh"
#include "G4ParallelNavigator.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelWorld::G4ParallelWorld(G4VPhysicalVolume &worldvolume)
 : fPstepper(new G4ParallelStepper),
   fPdriver(new G4ParallelNavigator(worldvolume))
{
  if (!fPstepper) {
    G4Exception("G4ParallelWorld::G4ParallelWorld()", "FatalException",
                FatalException, "NULL stepper. Failure in creation.");
  }
  if (!fPdriver) {
    G4Exception("G4ParallelWorld::G4ParallelWorld()", "FatalException",
                FatalException, "NULL navigator driver. Failure in creation.");
  }
}

G4ParallelWorld::~G4ParallelWorld()
{
  delete fPdriver;
  delete fPstepper;
}

G4VParallelStepper &G4ParallelWorld::GetParallelStepper()
{
  return *fPstepper;
}

G4VPGeoDriver &G4ParallelWorld::GetGeoDriver()
{
  return *fPdriver;
}
