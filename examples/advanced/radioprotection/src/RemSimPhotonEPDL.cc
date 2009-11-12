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
// $Id: RemSimPhotonEPDL.cc,v 1.8 2009-11-12 05:12:18 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#include "RemSimPhotonEPDL.hh"

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

#include "G4StepLimiter.hh"

RemSimPhotonEPDL::RemSimPhotonEPDL(const G4String& name):
 G4VPhysicsConstructor(name)
{ }

RemSimPhotonEPDL::~RemSimPhotonEPDL()
{ }

void RemSimPhotonEPDL::ConstructProcess()
{
  // Add EPDL processes for photons
  
  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* pmanager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
     
      if (particleName == "gamma") 
	{
	  // Photon     
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      theRayleigh->SetModel(new G4LivermoreRayleighModel());  //not strictly necessary
      pmanager->AddDiscreteProcess(theRayleigh);

      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetModel(new G4LivermorePhotoElectricModel());
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
	
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetModel(new G4LivermoreComptonModel());
      pmanager->AddDiscreteProcess(theComptonScattering);
	
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetModel(new G4LivermoreGammaConversionModel());
      pmanager->AddDiscreteProcess(theGammaConversion);  
     
      pmanager->AddProcess(new G4StepLimiter(),  -1,-1,3);
	}   
    }
}
