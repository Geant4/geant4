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
//
// LocalPenelopeEmPhysic.cc
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
//
// See more at: http://workgroup.lngs.infn.it/geant4lns/
// ----------------------------------------------------------------------------

#include "LocalPenelopeEmPhysic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4PenelopeCompton.hh"
#include "G4PenelopeGammaConversion.hh"
#include "G4PenelopePhotoElectric.hh"
#include "G4PenelopeRayleigh.hh"

#include "G4MultipleScattering.hh"

#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LocalPenelopeEmPhysic::LocalPenelopeEmPhysic(const G4String& name)
  :  G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LocalPenelopeEmPhysic::~LocalPenelopeEmPhysic()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LocalPenelopeEmPhysic::ConstructProcess()
{
  // Add Penelope or standard EM Processes

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma         
      pmanager->AddDiscreteProcess(new G4PenelopePhotoElectric);
      pmanager->AddDiscreteProcess(new G4PenelopeCompton);
      pmanager->AddDiscreteProcess(new G4PenelopeGammaConversion);
      pmanager->AddDiscreteProcess(new G4PenelopeRayleigh);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering,     -1, 1, 1);
      pmanager->AddProcess(new G4PenelopeIonisation,     -1, 2, 2);
      pmanager->AddProcess(new G4PenelopeBremsstrahlung, -1,-1, 3);
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering,     -1, 1, 1);
      pmanager->AddProcess(new G4PenelopeIonisation,     -1, 2, 2);
      pmanager->AddProcess(new G4PenelopeBremsstrahlung, -1,-1, 3);      
      pmanager->AddProcess(new G4PenelopeAnnihilation,    0,-1, 4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 4, 4);       
     
    } else if( particleName == "alpha" || particleName == "GenericIon" ) { 
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,        -1, 2, 2);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

