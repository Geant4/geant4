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
//---------------------------------------------------------------------------
//
// ClassName:   G4QAtomicPhysics
//
// Author:      M. Kosov 20.11.2009 (similar to G4QAtomicPhysics)
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QAtomicPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiOmegaMinus.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4BuilderType.hh"

G4QAtomicPhysics::G4QAtomicPhysics(G4int ver)
  : G4VPhysicsConstructor("CHIPS Atomic"), verbose(ver)
{
  G4LossTableManager::Instance();
  SetPhysicsType(bElectromagnetic);
}

G4QAtomicPhysics::G4QAtomicPhysics(G4int ver, const G4String& name)
  : G4VPhysicsConstructor(name), verbose(ver)
{
  G4LossTableManager::Instance();
  SetPhysicsType(bElectromagnetic);
}

G4QAtomicPhysics::~G4QAtomicPhysics()
{}

void G4QAtomicPhysics::ConstructParticle()
{
// gamma
  G4Gamma::Gamma();

// leptons
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();

// mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();

// barions
  G4Proton::Proton();
  G4AntiProton::AntiProton();
  G4SigmaPlus::SigmaPlus();
  G4SigmaMinus::SigmaMinus();
  G4XiMinus::XiMinus();
  G4OmegaMinus::OmegaMinus();
  G4AntiSigmaPlus::AntiSigmaPlus();
  G4AntiSigmaMinus::AntiSigmaMinus();
  G4AntiXiMinus::AntiXiMinus();
  G4AntiOmegaMinus::AntiOmegaMinus();

// ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIonDefinition();
}

void G4QAtomicPhysics::ConstructProcess()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particlePointer = theParticleIterator->value();
    G4ProcessManager* pmanager = particlePointer->GetProcessManager();
    G4String particle = particlePointer->GetParticleName();
    if(verbose > 1) G4cout<<"###G4QAtomicPhysics::ConstructProcesses: try "
                          <<GetPhysicsName()<<" builder for "<<particle<<G4endl;
    if     ( particle == "gamma")
    {

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
    }
    else if( particle == "e-")
    {
      pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,           -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),     -1,-3, 3);
    }
    else if( particle == "e+")
    {
      pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,           -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,       -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,      0,-1, 4);
    }
    else if( particle == "mu+"      || particle == "mu-" )
    {
      pmanager->AddProcess(new G4MuMultipleScattering,  -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,          -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,      -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,      -1,-4, 4);
    }
    else if( particle == "pi-"      || particle == "pi+"         ||
             particle == "kaon-"    || particle == "kaon+"       ||
             particle == "proton"   || particle == "anti_proton" ||
             particle == "tau-"     || particle == "tau+"        ||
             particle == "deuteron" || particle == "triton"      ||
             particle == "xi-"      || particle == "anti_xi-"    ||
             particle == "sigma+"   || particle == "anti_sigma+" ||
             particle == "sigma-"   || particle == "anti_sigma-" ||
             particle == "omega-"   || particle == "anti_omega-" )
    {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
    }
    else if(particle == "alpha" || particle == "He3" ||particle == "GenericIon")
    {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
    }
    else if (particle == "B+"        || particle == "B-"             ||
             particle == "D+"        || particle == "D-"             ||
             particle == "Ds+"       || particle == "Ds-"            ||
             particle == "lambda_c+" || particle == "anti_lambda_c+" ||
             particle == "sigma_c+"  || particle == "anti_sigma_c+"  ||
             particle == "sigma_c++" || particle == "anti_sigma_c++" ||
             particle == "xi_c+"     || particle == "anti_xi_c+"     )
    {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
    }
  }
  G4EmProcessOptions opt;
  opt.SetVerbose(verbose);
}
