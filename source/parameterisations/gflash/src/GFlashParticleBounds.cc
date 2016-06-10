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
// $Id: GFlashParticleBounds.cc 68057 2013-03-13 14:46:00Z gcosmo $
//
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ---------------- GFlashParticleBounds ----------------
//
// Author: Joanna Weng - 9.11.2004
// ------------------------------------------------------------

#include "GFlashParticleBounds.hh"

#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


GFlashParticleBounds::GFlashParticleBounds()
{    
  // e+e- defaults
  EMinEneToParametrise = 0.10*GeV;  
  EMaxEneToParametrise = 10000.00*GeV;  
  EEneToKill = 0.1*GeV; // Energie at which electrons are killed
}

GFlashParticleBounds::~GFlashParticleBounds()
{
}

void GFlashParticleBounds::
SetMinEneToParametrise(G4ParticleDefinition &particleType, G4double enemin)
{ 
  if( &particleType == G4Electron::ElectronDefinition()||
    &particleType == G4Positron::PositronDefinition()) 
  EMinEneToParametrise = enemin;
}

void GFlashParticleBounds::
SetMaxEneToParametrise(G4ParticleDefinition &particleType, G4double enemax)
{
  if( &particleType == G4Electron::ElectronDefinition()||
    &particleType == G4Positron::PositronDefinition()) 
  EMaxEneToParametrise = enemax; 
}

void GFlashParticleBounds::
SetEneToKill(G4ParticleDefinition &particleType, G4double enekill)
{
  if( &particleType == G4Electron::ElectronDefinition()||
    &particleType == G4Positron::PositronDefinition()) 
  EEneToKill = enekill; 
}

G4double GFlashParticleBounds::
GetMinEneToParametrise(G4ParticleDefinition &particleType) 
{ 
  G4double result = DBL_MAX;
  if( &particleType == G4Electron::ElectronDefinition()||
    &particleType == G4Positron::PositronDefinition()) 
  {
    result = EMinEneToParametrise;
  }        
  return result;
}

G4double GFlashParticleBounds::
GetMaxEneToParametrise(G4ParticleDefinition &particleType) 
{ 
  G4double result = 0;
  if( &particleType == G4Electron::ElectronDefinition()||
    &particleType == G4Positron::PositronDefinition()) 
  {
    result = EMaxEneToParametrise; 
  }
  return result;
}

G4double GFlashParticleBounds::
GetEneToKill(G4ParticleDefinition & particleType) 
{
  if (&particleType == G4Electron::ElectronDefinition() ||
    &particleType == G4Positron::PositronDefinition())  
  return EEneToKill;
  else return (-DBL_MAX);
}
