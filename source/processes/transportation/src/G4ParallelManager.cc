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
// $Id: G4ParallelManager.cc,v 1.2 2002-04-09 17:40:16 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelManager.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelManager.hh"
#include "G4ParallelStepper.hh"
#include "G4ParallelWorld.hh"
#include "G4ParallelNavigator.hh"
#include "G4ParallelTransport.hh"
#include "G4ProcessPlacer.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelManager::G4ParallelManager(G4VPhysicalVolume &worldvolume,
				     const G4String &particlename)
 : fPworld(new G4ParallelWorld(worldvolume)),
   fParticleName(particlename),
   fParallelTransport(0)
{}

  
G4ParallelManager::~G4ParallelManager()
{
  delete fPworld;
  if (fParallelTransport) delete fParallelTransport;
}

G4ParallelWorld &G4ParallelManager::GetParallelWorld()
{
  return *fPworld;
}

G4String G4ParallelManager::GetParticleName()
{
  return fParticleName;
}

G4ParallelTransport *G4ParallelManager::CreateParallelTransport()
{
  if (!fParallelTransport) {
    fParallelTransport = new G4ParallelTransport(fPworld->GetGeoDriver(), 
				       fPworld->GetParallelStepper());
  }
  return fParallelTransport;
}

void G4ParallelManager::Initialize()
{
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateParallelTransport());
}
