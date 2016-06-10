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
// $Id: G4NeutronKiller.cc 68048 2013-03-13 14:34:07Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronKiller
//
// Description: The process to kill particles to save CPU
//
// Author:      V.Ivanchenko 26/09/00 for HARP software
//
//----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4NeutronKiller.hh"

#include "G4SystemOfUnits.hh"
#include "G4NeutronKillerMessenger.hh"
#include "G4TransportationProcessType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NeutronKiller::G4NeutronKiller(const G4String& processName, G4ProcessType aType)
 : G4VDiscreteProcess(processName, aType)
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(NEUTRON_KILLER));

  kinEnergyThreshold = 0.0;
  timeThreshold = DBL_MAX;
  pMess = new G4NeutronKillerMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NeutronKiller::~G4NeutronKiller() 
{ 
  delete pMess; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4NeutronKiller::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetParticleName() == "neutron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NeutronKiller::SetTimeLimit(G4double val)
{
  timeThreshold = val;
  if(verboseLevel > 0) 
    G4cout << "### G4NeutronKiller: timeLimit(ns) = " 
	   << timeThreshold/ns << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NeutronKiller::SetKinEnergyLimit(G4double val)
{
  kinEnergyThreshold = val;
  if(verboseLevel > 0) 
    G4cout << "### G4NeutronKiller: Tracking cut E(MeV) = " 
	   << kinEnergyThreshold/MeV << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


