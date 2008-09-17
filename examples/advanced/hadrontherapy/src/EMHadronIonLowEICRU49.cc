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
// $Id: EMHadronIonLowEICRU49.cc; 
// Last modified: A.Lechner (anton.lechner@cern.ch), August 2008;
//
// See more at: http://geant4infn.wikispaces.com
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------


#include "G4SDManager.hh"
#include "EMHadronIonLowEICRU49.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4StepLimiter.hh"


EMHadronIonLowEICRU49::EMHadronIonLowEICRU49(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (ions and charged hadrons)" 
        << G4endl
        << "                             G4hLowEnergyIonisation (ions and charged hadrons)" 
        << G4endl
        << "                             (applying ICRU49 stopping power parametrisations)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMHadronIonLowEICRU49::~EMHadronIonLowEICRU49()
{ }

void EMHadronIonLowEICRU49::ConstructProcess()
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
      
      if (charge != 0.0         && 
          particleName != "e+"  && 
          particleName != "mu+" &&
	  particleName != "e-"  && 
          particleName != "mu-") 
	{
	  if((!particle -> IsShortLived()) &&
	     (particle -> GetParticleName() != "chargedgeantino"))
	    {	  

              G4MultipleScattering* hadronIonMultipleScatProcess = new G4MultipleScattering(); 
              G4hLowEnergyIonisation* hadronIonIonisationProcess = new G4hLowEnergyIonisation();
              hadronIonIonisationProcess -> SetNuclearStoppingOn();
	  
              G4StepLimiter* hadronIonStepLimiter = new G4StepLimiter();

              processManager -> AddProcess(hadronIonMultipleScatProcess, -1, 1, 1);   
	      processManager -> AddProcess(hadronIonIonisationProcess, -1, 2, 2);
              processManager -> AddProcess(hadronIonStepLimiter, -1, -1, 3);
	    }
	}
    }
}

   
