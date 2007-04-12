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
// $Id: Tst52PhotonEPDL.cc,v 1.2 2007-04-12 14:51:52 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: guatelli@ge.infn.it, Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst52PhotonEPDL.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4StepLimiter.hh"

Tst52PhotonEPDL::Tst52PhotonEPDL(const G4String& name): G4VPhysicsConstructor(name)
{ 
  G4cout <<"Photon EPDL initialized!!!" << G4endl;
  cut = 0.001 * keV;
}

Tst52PhotonEPDL::~Tst52PhotonEPDL()
{ }

void Tst52PhotonEPDL::ConstructProcess()
{
  // Add EPDL processes for photons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "gamma") 
	{
          photo = new G4LowEnergyPhotoElectric();
	  manager->AddDiscreteProcess(photo);
	  manager->AddDiscreteProcess(new G4LowEnergyCompton);
	  manager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
	  manager->AddDiscreteProcess(new G4LowEnergyRayleigh);
	  manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}   
    }
}
