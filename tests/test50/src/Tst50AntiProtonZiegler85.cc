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
// $Id: Tst50AntiProtonZiegler85.cc,v 1.3 2010-06-25 09:46:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Stephane Chauvie (chauvie@to.infn.it)
//
// History:
// -----------
// 22 Feb 2006         SC      Created
//
// -------------------------------------------------------------------

#include "Tst50AntiProtonZiegler85.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
//#include "G4MultipleScattering.hh"
#include "G4AntiProton.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hLowEnergyLoss.hh"
#include "G4hZiegler1985p.hh"
#include "G4StepLimiter.hh"

Tst50AntiProtonZiegler85::Tst50AntiProtonZiegler85(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50AntiProtonZiegler85::~Tst50AntiProtonZiegler85()
{ }

void Tst50AntiProtonZiegler85::ConstructProcess()
{

  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "anti_proton") 
	{
	  // G4VProcess*  multipleScattering= new G4MultipleScattering(); 
	  G4hLowEnergyIonisation* ion= new G4hLowEnergyIonisation();
	  
	  ion -> SetElectronicStoppingPowerModel(particle, "Ziegler1985p");
          ion -> SetNuclearStoppingPowerModel("Ziegler1985");
          ion -> SetNuclearStoppingOn() ;
          ion -> SetBarkasOn() ;
	  //  manager->AddProcess(multipleScattering,-1,1,1);  	
	  ion -> SetEnlossFluc(false);          
          manager -> AddProcess(ion,-1,2,2);
          manager -> AddProcess (new G4StepLimiter(),-1,-1,3);
	}
    }
}
