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
// $Id: G4ParallelWorld.cc,v 1.2 2002-04-09 16:23:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
 : fWorldVolume(&worldvolume),
   fPstepper(new G4ParallelStepper),
   fPdriver(new G4ParallelNavigator(*fWorldVolume))
{}

G4ParallelWorld::~G4ParallelWorld()
{
  delete fPdriver;
  delete fPstepper;
}

G4ParallelWorld::G4ParallelWorld(const G4ParallelWorld &rhs)
 : fWorldVolume(rhs.GetWorldVolume())
{
   fPstepper = new G4ParallelStepper;
   fPdriver = new G4ParallelNavigator(*fWorldVolume);
}

G4ParallelWorld &G4ParallelWorld::operator=(const G4ParallelWorld &rhs)
{
  if (this!=&rhs) {

    fWorldVolume = rhs.GetWorldVolume();
    fPstepper = new G4ParallelStepper;
    fPdriver = new G4ParallelNavigator(*fWorldVolume);
  }
  return *this;
}


G4VPhysicalVolume *G4ParallelWorld::GetWorldVolume() const
{
  return fWorldVolume;
}


G4VParallelStepper &G4ParallelWorld::GetParallelStepper()
{
  return *fPstepper;
}

G4VPGeoDriver &G4ParallelWorld::GetGeoDriver()
{
  return *fPdriver;
}
