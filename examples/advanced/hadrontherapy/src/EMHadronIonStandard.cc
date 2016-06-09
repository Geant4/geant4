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
// $Id: EMIonStandard.cc; November 2008
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
// 
// This class manages the electromagnetic processes for charged hadrons and ions
// using the Standard Models of Geant4
// ----------------------------------------------------------------------------

#include "EMHadronIonStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4MultipleScattering.hh"
#include "G4StepLimiter.hh"
#include "G4EmProcessOptions.hh"
#include "G4MuIonisation.hh"

EMHadronIonStandard::EMHadronIonStandard(const G4String& name): 
   G4VPhysicsConstructor(name)
{ }

EMHadronIonStandard::~EMHadronIonStandard()
{ }

void EMHadronIonStandard::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;

  // ********************************
  // *** Charged Hadrons and Ions ***
  // ********************************

  theParticleIterator -> reset();
  
  while( (*theParticleIterator)() )
    {
      particle = theParticleIterator -> value();
      processManager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
      G4double charge = particle -> GetPDGCharge();
  
      if (particleName == "GenericIon"|| 
	  particleName == "alpha" ||
	  particleName == "He3")
	{
	  processManager -> AddProcess(new G4hMultipleScattering, -1, 1, 1);   
	  processManager -> AddProcess(new G4hIonisation ,-1, 2, 2);
	}
      else
	{
	  if (charge != 0.0         && 
              particleName != "e+"  && 
              particleName != "mu+" &&
	      particleName != "e-"  && 
              particleName != "mu-") 
	    {
	      if((!particle -> IsShortLived()) &&
		 (particle -> GetParticleName() != "chargedgeantino"))
		{
		  processManager -> AddProcess(new G4hMultipleScattering, -1, 1, 1);   
		  processManager -> AddProcess(new G4hIonisation, -1, 2, 2);
		}
	    }
	}
    }

 }

   
