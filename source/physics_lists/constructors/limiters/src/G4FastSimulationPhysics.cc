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
//---------------------------------------------------------------------------
//
// ClassName:   G4FastSimulationPhysics
//
// Author:      M. Verderi (Nov.03.2016)
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FastSimulationPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4FastSimulationHelper.hh"
#include "G4FastSimulationManagerProcess.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4FastSimulationPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FastSimulationPhysics::G4FastSimulationPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name),
     fVerbose(false)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FastSimulationPhysics::~G4FastSimulationPhysics()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4FastSimulationPhysics::ActivateFastSimulation(const G4String particleName, const G4String parallelGeometryName)
{
  fParticlesUnderFastSimulation.push_back(particleName);
  fGeometries                  .push_back(parallelGeometryName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4FastSimulationPhysics::ConstructParticle()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4FastSimulationPhysics::ConstructProcess()
{
  
  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();
  
  while ( (*myParticleIterator)() )
    {
      G4ParticleDefinition*     particle = myParticleIterator->value();
      G4String              particleName = particle->GetParticleName();
      G4ProcessManager*         pmanager = particle->GetProcessManager();
      
      // -- include fast simulation manager process interface:
      auto itr = std::find( fParticlesUnderFastSimulation.begin(),
			    fParticlesUnderFastSimulation.end(),
			    particleName                           );
      
      if ( itr != fParticlesUnderFastSimulation.end() )
	{
	  std::size_t ipos = itr - fParticlesUnderFastSimulation.begin();
	  G4String geometry = fGeometries[ipos];
	  if ( geometry == "" ) G4FastSimulationHelper::ActivateFastSimulation(pmanager);
	  else                  G4FastSimulationHelper::ActivateFastSimulation(pmanager, geometry);
	}
    }

  // -- tells what is done:
  if ( fVerbose )
    {
      // -- print:
      myParticleIterator->reset();
      
      while ( (*myParticleIterator)() )
	{
	  G4ParticleDefinition*     particle = myParticleIterator->value();
	  G4String              particleName = particle->GetParticleName();
	  G4ProcessManager*         pmanager = particle->GetProcessManager();
	  
	  G4bool isUnderFastSimulation(false);
	  G4String processAndGeometryNames;
	  G4int icount(0);
	  
	  G4ProcessVector*  vprocess = pmanager->GetProcessList();
	  for (G4int ip = 0 ; ip < (G4int)vprocess->size() ; ++ip)
	    {
	      G4VProcess* process = (*vprocess)[ip];
	      G4FastSimulationManagerProcess* pb = dynamic_cast< G4FastSimulationManagerProcess* >(process);
	      if ( pb != nullptr )
		{
		  isUnderFastSimulation = true;
		  if ( icount < 3 )
		    {
		      processAndGeometryNames += pb->GetProcessName();
		      processAndGeometryNames += "[geom:";
		      processAndGeometryNames += pb->GetWorldVolume()->GetName();
		      processAndGeometryNames += "] ";
		    }
		  else
		    {
		      processAndGeometryNames += "\n                 ";
		      processAndGeometryNames += pb->GetProcessName();
		      processAndGeometryNames += "[geom:";
		      processAndGeometryNames += pb->GetWorldVolume()->GetName();
		      processAndGeometryNames += "] ";
		      icount = 0;
		    }
		}
	    }
	  if ( isUnderFastSimulation ) G4cout << std::setw(14) << particleName << " : " << processAndGeometryNames << G4endl;
	}
    }
}
