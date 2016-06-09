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
// $Id: ExN07SteppingVerbose.cc,v 1.1 2006-11-04 19:23:07 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN07SteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4RunManagerKernel.hh"
#include "G4TrackingManager.hh"

ExN07SteppingVerbose::ExN07SteppingVerbose()
: nTimers(0),regIdx(-1),ep(false)
{}

ExN07SteppingVerbose::~ExN07SteppingVerbose()
{
  for(G4int j=0;j<nTimers;j++)
  { delete fTimers[j]; }
  fTimers.clear();
}

void ExN07SteppingVerbose::InitializeTimers()
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  nRegions = regionStore->size();
  nTimers = 2 * nRegions;
  G4int nEnt = fTimers.size();
  if(nEnt<nTimers)
  {
    for(G4int i=nEnt;i<nTimers;i++)
    { fTimers.push_back(new G4SliceTimer); }
  }
  for(G4int j=0;j<nTimers;j++)
  { fTimers[j]->Clear(); }
  regIdx = -1;
  ep = false;

  // Set verbosity for timing
  G4RunManagerKernel::GetRunManagerKernel()->GetTrackingManager()->SetVerboseLevel(0);
  fManager->SetVerboseLevel(1);
}

void ExN07SteppingVerbose::Report()
{
  for(G4int i=0;i<nRegions;i++)
  {
    G4cout << G4endl;
    G4cout << "Region <" << (*G4RegionStore::GetInstance())[i]->GetName() << ">" << G4endl;
    G4cout << " All particles : User=" << fTimers[i]->GetUserElapsed()
     << "  Real=" << fTimers[i]->GetRealElapsed()
     << "  Sys=" << fTimers[i]->GetSystemElapsed() << G4endl;
    G4cout << " e+ / e-       : User=" << fTimers[nRegions+i]->GetUserElapsed()
     << "  Real=" << fTimers[nRegions+i]->GetRealElapsed()
     << "  Sys=" << fTimers[nRegions+i]->GetSystemElapsed() << G4endl;
  }
  G4cout << G4endl;
}

void ExN07SteppingVerbose::NewStep()
{
  CopyState();
  G4Region* reg = fTrack->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion();
  regIdx = FindRegion(reg);
  fTimers[regIdx]->Start();
  G4ParticleDefinition* pd = fTrack->GetDefinition();
  if(pd==G4Electron::ElectronDefinition() || pd==G4Positron::PositronDefinition()) ep = true;
  if(ep) fTimers[nRegions+regIdx]->Start();
} 

void ExN07SteppingVerbose::StepInfo()
{
  fTimers[regIdx]->Stop();
  if(ep)
  {
    fTimers[nRegions+regIdx]->Stop();
    ep = false;
  }
  regIdx = -1;
}

G4int ExN07SteppingVerbose::FindRegion(G4Region* rgn)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4int sz = regionStore->size();
  for(G4int i=0;i<sz;i++)
  { if(rgn==(*regionStore)[i]) return i; }
  return -1;
}


// Empty methods not to be used.
void ExN07SteppingVerbose::TrackBanner() 
{;}
void ExN07SteppingVerbose::AtRestDoItInvoked() 
{;}
void ExN07SteppingVerbose::AlongStepDoItAllDone()
{;}
void ExN07SteppingVerbose::PostStepDoItAllDone()
{;}
void ExN07SteppingVerbose::AlongStepDoItOneByOne()
{;}
void ExN07SteppingVerbose::PostStepDoItOneByOne()
{;}
void ExN07SteppingVerbose::TrackingStarted()
{;}
void ExN07SteppingVerbose::DPSLStarted()
{;}
void ExN07SteppingVerbose::DPSLUserLimit()
{;}
void ExN07SteppingVerbose::DPSLPostStep()
{;}
void ExN07SteppingVerbose::DPSLAlongStep()
{;}
void ExN07SteppingVerbose::VerboseTrack()
{;}
void ExN07SteppingVerbose::VerboseParticleChange()
{;}


