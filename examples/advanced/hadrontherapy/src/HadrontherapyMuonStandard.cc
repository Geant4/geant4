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
// $Id: HadrontherapyMuonStandard.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyMuonStandard.hh"
#include "G4ProcessManager.hh"

//muon:
#include "G4hIonisation.hh"
#include "G4MultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh" 
#include "G4StepLimiter.hh"

HadrontherapyMuonStandard::HadrontherapyMuonStandard(const G4String& name): G4VPhysicsConstructor(name)
{ }

HadrontherapyMuonStandard::~HadrontherapyMuonStandard()
{ }

void HadrontherapyMuonStandard::ConstructProcess()
{
  // Add standard processes for muons
  
  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
     
      if(( particleName == "mu+")|| (particleName == "mu-" ))
	{//muon
	  G4VProcess* aMultipleScattering = new G4MultipleScattering();
	  G4VProcess* aBremsstrahlung     = new G4MuBremsstrahlung();
	  G4VProcess* aPairProduction     = new G4MuPairProduction();
	  G4VProcess* anIonisation        = new G4MuIonisation();
	  //
	  // add processes
	  pmanager -> AddProcess(anIonisation);
	  pmanager -> AddProcess(aMultipleScattering);
	  pmanager -> AddProcess(aBremsstrahlung);
	  pmanager -> AddProcess(aPairProduction);
	  //
	  // set ordering for AlongStepDoIt
	  pmanager -> SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
	  pmanager -> SetProcessOrdering(anIonisation,        idxAlongStep,2);
	  pmanager -> SetProcessOrdering(aBremsstrahlung,     idxAlongStep,3);
	  pmanager -> SetProcessOrdering(aPairProduction,     idxAlongStep,4);      
	  //
	  // set ordering for PostStepDoIt
	  pmanager -> SetProcessOrdering(aMultipleScattering, idxPostStep,1);
	  pmanager -> SetProcessOrdering(anIonisation,        idxPostStep,2);
	  pmanager -> SetProcessOrdering(aBremsstrahlung,     idxPostStep,3);
	  pmanager -> SetProcessOrdering(aPairProduction,     idxPostStep,4);
          pmanager -> AddProcess(new G4StepLimiter(),      -1,-1, 3);

	  if( particleName == "mu-" )
	    pmanager -> AddProcess(new G4MuonMinusCaptureAtRest(), 0,-1,-1);
	}
    }
}
