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
/// \file RE06/src/RE06SteppingVerbose.cc
/// \brief Implementation of the RE06SteppingVerbose class
//
//

#include "RE06SteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4RunManagerKernel.hh"
#include "G4TrackingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE06SteppingVerbose::RE06SteppingVerbose()
: G4VSteppingVerbose(),
  fTimers(),
  fNofTimers(0),
  fRegIdx(-1),
  fEp(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE06SteppingVerbose::~RE06SteppingVerbose()
{
  for(G4int j=0;j<fNofTimers;j++)
  { delete fTimers[j]; }
  fTimers.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06SteppingVerbose::InitializeTimers()
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  fNofRegions = regionStore->size();
  fNofTimers = 2 * fNofRegions;
  G4int nEnt = fTimers.size();
  if(nEnt<fNofTimers)
  {
    for(G4int i=nEnt;i<fNofTimers;i++)
    { fTimers.push_back(new G4SliceTimer); }
  }
  for(G4int j=0;j<fNofTimers;j++)
  { fTimers[j]->Clear(); }
  fRegIdx = -1;
  fEp = false;

  // Set verbosity for timing
  G4RunManagerKernel::GetRunManagerKernel()->GetTrackingManager()->SetVerboseLevel(0);
#ifdef G4VERBOSE
  fManager->SetVerboseLevel(1);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06SteppingVerbose::Report()
{
  for(G4int i=0;i<fNofRegions;i++)
  {
    G4cout << G4endl;
    G4cout << "Region <" 
     << (*G4RegionStore::GetInstance())[i]->GetName() << ">" << G4endl;
    G4cout << " All particles : User=" << fTimers[i]->GetUserElapsed()
     << "  Real=" << fTimers[i]->GetRealElapsed()
     << "  Sys=" << fTimers[i]->GetSystemElapsed() << G4endl;
    G4cout << " e+ / e-       : User=" << fTimers[fNofRegions+i]->GetUserElapsed()
     << "  Real=" << fTimers[fNofRegions+i]->GetRealElapsed()
     << "  Sys=" << fTimers[fNofRegions+i]->GetSystemElapsed() << G4endl;
  }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06SteppingVerbose::NewStep()
{
  CopyState();
  G4Region* reg = fTrack->GetStep()->GetPreStepPoint()
                  ->GetPhysicalVolume()->GetLogicalVolume()->GetRegion();
  fRegIdx = FindRegion(reg);
  fTimers[fRegIdx]->Start();
  G4ParticleDefinition* pd = fTrack->GetDefinition();
  if(pd==G4Electron::ElectronDefinition() || 
     pd==G4Positron::PositronDefinition()) fEp = true;
  if(fEp) fTimers[fNofRegions+fRegIdx]->Start();
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06SteppingVerbose::StepInfo()
{
  fTimers[fRegIdx]->Stop();
  if(fEp)
  {
    fTimers[fNofRegions+fRegIdx]->Stop();
    fEp = false;
  }
  fRegIdx = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int RE06SteppingVerbose::FindRegion(G4Region* rgn)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4int sz = regionStore->size();
  for(G4int i=0;i<sz;i++)
  { if(rgn==(*regionStore)[i]) return i; }
  return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
