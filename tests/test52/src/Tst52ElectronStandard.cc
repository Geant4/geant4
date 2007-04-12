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
// $Id: Tst52ElectronStandard.cc,v 1.1 2007-04-12 12:00:17 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:guatelli@ge.infn.it, Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst52ElectronStandard.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4StepLimiter.hh"

Tst52ElectronStandard::Tst52ElectronStandard(const G4String& name): G4VPhysicsConstructor(name)
{ G4cout<< "Electron Standard Processes initialized !!!!" << G4endl;
 facValue = 0.02;//default value
}

Tst52ElectronStandard::~Tst52ElectronStandard()
{ }

void Tst52ElectronStandard::ConstructProcess()
{
  // Add standard processes for electrons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "e-") 
	{
          scattering = new G4MultipleScattering();
	  manager->AddProcess(scattering, -1, 1,1);
	  manager -> AddProcess(new G4eIonisation,        -1, 2,2);
	  manager -> AddProcess(new G4eBremsstrahlung,    -1,-1,3);
     
	  // max step
	  manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	  
	  // Facrange
	  scattering -> SetFacrange(facValue);
	}   
    }
}
void Tst52ElectronStandard::SetFacRange(G4double value)
{
 facValue = value;
 G4cout << "The value of the facRange is:" << value << " for electrons"<< G4endl;
}
