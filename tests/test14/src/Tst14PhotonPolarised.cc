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
// $Id: Tst14PhotonPolarised.cc,v 1.5 2010-08-28 20:35:36 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst14PhotonPolarised.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
/*
#include "G4LowEnergyPolarizedCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyPolarizedRayleigh.hh"
*/


#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePolarizedPhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermorePolarizedRayleighModel.hh"


Tst14PhotonPolarised::Tst14PhotonPolarised(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst14PhotonPolarised::~Tst14PhotonPolarised()
{ }

void Tst14PhotonPolarised::ConstructProcess()
{
  // Add polarised processes for photons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "gamma") 
	{
/*
	  manager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
	  manager->AddDiscreteProcess(new G4LowEnergyPolarizedCompton);
	  manager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
	  manager->AddDiscreteProcess(new G4LowEnergyPolarizedRayleigh);
*/


      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePolarizedPhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePolarizedPhotoElectricModel();
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      manager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermorePolarizedComptonModel* theLivermoreComptonModel = new G4LivermorePolarizedComptonModel();
      theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
      manager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermorePolarizedGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermorePolarizedGammaConversionModel();
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      manager->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermorePolarizedRayleighModel* theRayleighModel = new G4LivermorePolarizedRayleighModel();
      theRayleigh->AddEmModel(0, theRayleighModel);
      manager->AddDiscreteProcess(theRayleigh);



	}   
    }
}
