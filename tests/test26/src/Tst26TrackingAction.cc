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
// $Id: Tst26TrackingAction.cc,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26TrackingAction.hh"
#include "Tst26RunAction.hh"

#include "G4Track.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26TrackingAction::Tst26TrackingAction(Tst26RunAction* run)
  :Tst26Run(run),
   vertex(0),
   muon(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  if(!vertex) vertex = (G4RegionStore::GetInstance())->GetRegion("VertexDetector");
  if(!muon)   muon = (G4RegionStore::GetInstance())->GetRegion("MuonDetector");
  if(1 == track->GetTrackID()) return;
  G4int regionIndex = 0;
  const G4Region* r = track->GetVolume()->GetLogicalVolume()->GetRegion();
  if(r == vertex) regionIndex = 1;
  else if(r == muon) regionIndex = 2;
  G4int particleIndex = -1;
  const G4ParticleDefinition* pd = track->GetDefinition();
  if(G4Gamma::Gamma() == pd) particleIndex = 0;
  else if(G4Electron::Electron() == pd) particleIndex = 1;
  else if(G4Positron::Positron() == pd) particleIndex = 2;
 
  if(particleIndex >= 0) Tst26Run->AddParticle(particleIndex,regionIndex); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


