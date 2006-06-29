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
// $Id: Tst26TrackingAction.cc,v 1.6 2006-06-29 21:54:05 gunter Exp $
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


