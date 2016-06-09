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
// $Id: RemSimPrimaryGeneratorAction.cc,v 1.14 2005/12/02 10:08:34 guatelli Exp $// Author: Susanna Guatelli, guatelli@ge.infn.it

#include "RemSimPrimaryGeneratorAction.hh"
#ifdef G4ANALYSIS_USE  
#include "RemSimAnalysisManager.hh"
#endif
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
RemSimPrimaryGeneratorAction::RemSimPrimaryGeneratorAction()
{
  particleGun = new G4GeneralParticleSource();
}

RemSimPrimaryGeneratorAction::~RemSimPrimaryGeneratorAction()
{
   delete particleGun;
}

G4double RemSimPrimaryGeneratorAction::GetInitialEnergy()
{
   G4double primaryParticleEnergy = particleGun -> GetParticleEnergy();
   return primaryParticleEnergy;   
}

void RemSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 
#ifdef G4ANALYSIS_USE   
 G4double energy = particleGun -> GetParticleEnergy(); 
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> primaryParticleEnergyDistribution(energy/MeV);
#endif 

particleGun -> GeneratePrimaryVertex(anEvent);
}
