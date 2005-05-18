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
//    **************************************
//    *                                    *
//    *    HadrontherapyProtonHadro.cc        *
//    *                                    *
//    **************************************
//
// $Id: HadrontherapyProtonHadro.cc,v 1.3 2005-05-18 07:53:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 
#include "HadrontherapyProtonHadro.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4LElastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4ProtonInelasticProcess.hh"

HadrontherapyProtonHadro::HadrontherapyProtonHadro(const G4String& name): 
G4VPhysicsConstructor(name)
{
}
HadrontherapyProtonHadro::~HadrontherapyProtonHadro()
{}

void HadrontherapyProtonHadro::ConstructProcess()
{
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      G4LElastic* theElasticModel = new G4LElastic;
      theElasticProcess->RegisterMe(theElasticModel);
  
      theParticleIterator->reset();
      while ((*theParticleIterator)()) 
	{
	  G4ParticleDefinition* particle = theParticleIterator->value();
	  G4ProcessManager* pmanager = particle->GetProcessManager();
	  G4String particleName = particle->GetParticleName();
      
	  
	  if (particleName == "proton") 
	    {
	      pmanager->AddDiscreteProcess(theElasticProcess);
	      G4ProtonInelasticProcess* theInelasticProcess =
                                  new G4ProtonInelasticProcess("inelastic");
	      G4LEProtonInelastic* theLEInelasticModel = new G4LEProtonInelastic;
	      theInelasticProcess->RegisterMe(theLEInelasticModel);
	      G4HEProtonInelastic* theHEInelasticModel = new G4HEProtonInelastic;
	      theInelasticProcess->RegisterMe(theHEInelasticModel);
	      pmanager->AddDiscreteProcess(theInelasticProcess);
	    }
	}
}



