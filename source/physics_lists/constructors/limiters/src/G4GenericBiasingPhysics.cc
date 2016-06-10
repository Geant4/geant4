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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4GenericBiasingPhysics
//
// Author:      M. Verderi (Sept.10.2013)
// Modified:
// 07/11/2014, M. Verderi : fix bug of PhysicsBias(...) which was not taking
//             into account the vector of processes passed, but biasing all.
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GenericBiasingPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4BiasingHelper.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4GenericBiasingPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4GenericBiasingPhysics::G4GenericBiasingPhysics(const G4String& name)
   :  G4VPhysicsConstructor(name)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4GenericBiasingPhysics::~G4GenericBiasingPhysics()
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::PhysicsBias(const G4String& particleName)
{
  fBiasedParticles.push_back(particleName);
  std::vector< G4String > dummy;
  fBiasedProcesses.push_back(dummy);
  fBiasAllProcesses.push_back(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::PhysicsBias(const G4String& particleName, const std::vector< G4String >& processNames)
{
  fBiasedParticles.push_back(particleName);
  fBiasedProcesses.push_back(processNames);
  fBiasAllProcesses.push_back(false); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::NonPhysicsBias(const G4String& particleName)
{
  fNonPhysBiasedParticles.push_back(particleName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::Bias(const G4String& particleName)
{
  PhysicsBias(particleName);
  NonPhysicsBias(particleName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::Bias(const G4String& particleName, const std::vector< G4String >& processNames)
{
  PhysicsBias(particleName, processNames);
  NonPhysicsBias(particleName);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::ConstructParticle()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::ConstructProcess()
{
  aParticleIterator->reset();
  
  while( (*aParticleIterator)() )
    {
      G4ParticleDefinition*     particle = aParticleIterator->value();
      G4String              particleName = particle->GetParticleName();
      G4ProcessManager*         pmanager = particle->GetProcessManager();
      
      // -- include non physics process interface for biasing:
      if ( std::find(fNonPhysBiasedParticles.begin(),
		     fNonPhysBiasedParticles.end(),
		     particleName             )  != fNonPhysBiasedParticles.end() )
	{
	  G4BiasingHelper::ActivateNonPhysicsBiasing(pmanager);
	}
      
      // -- wrap biased physics processes, all processes or only user selected:
      std::vector< G4String >::const_iterator particleIt =
	std::find(fBiasedParticles.begin(),
		  fBiasedParticles.end(),
		  particleName             );
      if ( particleIt == fBiasedParticles.end() ) continue;
      
      std::vector < G4String >& biasedProcesses = fBiasedProcesses [ particleIt - fBiasedParticles.begin() ];
      G4bool biasAll                            = fBiasAllProcesses[ particleIt - fBiasedParticles.begin() ];
      
      if ( biasAll )
	{
	  G4ProcessVector*  vprocess = pmanager->GetProcessList();
	  for (G4int ip = 0 ; ip < vprocess->size() ; ip++)
	    {
	      G4VProcess* process = (*vprocess)[ip];
	      biasedProcesses.push_back( process->GetProcessName() );
	    }
	}
      
      G4bool restartLoop(true);
      while ( restartLoop )
	{
	  for (std::size_t ip = 0 ; ip < biasedProcesses.size() ; ip++)
	    {
	      G4bool activ = G4BiasingHelper::ActivatePhysicsBiasing(pmanager, biasedProcesses[ip] );
	      restartLoop = activ;
	      if ( restartLoop ) break;
	    }
	}
      
    }
}

