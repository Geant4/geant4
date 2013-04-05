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
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst14PhotonEPDL.hh"

#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"

Tst14PhotonEPDL::Tst14PhotonEPDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst14PhotonEPDL::~Tst14PhotonEPDL()
{ }

void Tst14PhotonEPDL::ConstructProcess()
{
  // Add EPDL processes for photons
  
  aParticleIterator->reset();

  while( (*aParticleIterator)() )
    {
      G4ParticleDefinition* particle = aParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      G4double LivermoreHighEnergyLimit = GeV;

      if (particleName == "gamma") 
	{
	  G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
	  G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
	    new G4LivermorePhotoElectricModel();
	  theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
	  thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
	  manager->AddDiscreteProcess(thePhotoElectricEffect);

	  G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
	  G4LivermoreComptonModel* theLivermoreComptonModel = 
	    new G4LivermoreComptonModel();
	  theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
	  theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
	  manager->AddDiscreteProcess(theComptonScattering);

	  G4GammaConversion* theGammaConversion = new G4GammaConversion();
	  G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = 
	    new G4LivermoreGammaConversionModel();
	  theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
	  theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
	  manager->AddDiscreteProcess(theGammaConversion);

	  G4RayleighScattering* theRayleigh = new G4RayleighScattering();
	  G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
	  theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
	  theRayleigh->AddEmModel(0, theRayleighModel);
	  manager->AddDiscreteProcess(theRayleigh);
	}   
    }
}
