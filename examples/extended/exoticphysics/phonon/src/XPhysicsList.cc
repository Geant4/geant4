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
/// \file exoticphysics/phonon/src/XPhysicsList.cc
/// \brief Implementation of the XPhysicsList class
//
// $Id$
//

#include "XPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"              

#include "XTPhononSlow.hh"
#include "XTPhononFast.hh"
#include "XLPhonon.hh"

#include "XPhononScatteringProcess.hh"
#include "XPhononReflectionProcess.hh"
#include "XPhononDownconversionProcess.hh"

#include "G4SystemOfUnits.hh"


XPhysicsList::XPhysicsList():  G4VUserPhysicsList()
{
  G4cout<<"\n\nXPhysicsList::constructor: running"<<endl;
  defaultCutValue = 100*mm;
  SetVerboseLevel(1);
  G4cout<<"\n\nXPhysicsList::constructor: ran"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhysicsList::~XPhysicsList()
{}

void XPhysicsList::ConstructParticle()
{

  XLPhonon::PhononDefinition();
  XTPhononFast::PhononDefinition();
  XTPhononSlow::PhononDefinition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



void XPhysicsList::ConstructProcess()
{

  AddTransportation();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "XTPhononSlow") {
      G4cout<<"Registering Slow tansverse processes..."<<G4endl;     
      pmanager->AddDiscreteProcess(new XPhononScatteringProcess());      
      pmanager->AddDiscreteProcess(new XPhononDownconversionProcess());         
      pmanager->AddDiscreteProcess(new XPhononReflectionProcess());
    }
    
    if (particleName == "XTPhononFast") {
      
      G4cout<<"Registering fast transverse processes..."<<G4endl;
      pmanager->AddDiscreteProcess(new XPhononScatteringProcess());      
      pmanager->AddDiscreteProcess(new XPhononDownconversionProcess());         
      pmanager->AddDiscreteProcess(new XPhononReflectionProcess());
    }
        
    if (particleName == "XLPhonon") {      
      G4cout<<"Registering longitudinal processes..."<<G4endl;
      pmanager->AddDiscreteProcess(new XPhononScatteringProcess());
      pmanager->AddDiscreteProcess(new XPhononDownconversionProcess());        
      pmanager->AddDiscreteProcess(new XPhononReflectionProcess());      
    } 
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicsList::SetCuts()
{
  // These values are used as the default production thresholds
  // for the world volume.
  SetCutsWithDefault();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



