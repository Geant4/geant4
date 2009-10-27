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
// $Id: G4EmStandardPhysicsLYmsc.cc,v 1.3 2009-10-27 16:24:32 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysicsLYmsc
//
// Author:    V.Grichine based on Em_opt2 builder of  V.Ivanchenko 
//
// Modified:
// 19.12.2009 V.Grichine LY msc model was substituted
//
//----------------------------------------------------------------------------
//

#include "G4EmStandardPhysicsLYmsc.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4LYMultipleScattering.hh"
#include "G4CoulombScattering.hh"

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
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsLYmsc::G4EmStandardPhysicsLYmsc(
    G4int ver, const G4String& name)
  : G4VPhysicsConstructor(name), verbose(ver)
{
  G4LossTableManager::Instance();
  G4cout<<"G4EmStandardPhysicsLYmsc::G4EmStandardPhysicsLYmsc"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsLYmsc::~G4EmStandardPhysicsLYmsc()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsLYmsc::ConstructParticle()
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

// ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsLYmsc::ConstructProcess()
{
  // Add standard EM Processes
  G4cout<<"G4EmStandardPhysicsLYmsc::ConstructProcess()"<<G4endl;

  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if(verbose > 1)
    {
      G4cout << "### " << GetPhysicsName() << " instantiates for " 
	     << particleName << G4endl;
    }
    if (particleName == "gamma") 
    {
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
    } 
    else if (particleName == "e-") 
    {
      pmanager->AddProcess(new G4LYMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1,-3, 3);
      pmanager->AddDiscreteProcess(new G4CoulombScattering());
    } 
    else if (particleName == "e+") 
    {
      pmanager->AddProcess(new G4LYMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);
      pmanager->AddDiscreteProcess(new G4CoulombScattering());
    } 
    else if (particleName == "mu+" ||
             particleName == "mu-"    ) 
    {
      pmanager->AddProcess(new G4LYMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,    -1,-4, 4);
      pmanager->AddDiscreteProcess(new G4CoulombScattering());
    } 
    else if (particleName == "alpha" ||
             particleName == "He3" ||
             particleName == "GenericIon"   ) 
    {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
    } 
    else if (particleName == "pi+" ||
             particleName == "pi-" ||
	     particleName == "kaon+" ||
             particleName == "kaon-" ||
             particleName == "proton"    ) 
    {
      pmanager->AddProcess(new G4LYMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
    } 
    else if (  particleName == "B+" ||
	       particleName == "B-" ||
	       particleName == "D+" ||
	       particleName == "D-" ||
	       particleName == "Ds+" ||
	       particleName == "Ds-" ||
               particleName == "anti_lambda_c+" ||
               particleName == "anti_omega-" ||
               particleName == "anti_proton" ||
               particleName == "anti_sigma_c+" ||
               particleName == "anti_sigma_c++" ||
               particleName == "anti_sigma+" ||
               particleName == "anti_sigma-" ||
               particleName == "anti_xi_c+" ||
               particleName == "anti_xi-" ||
               particleName == "deuteron" ||
	       particleName == "lambda_c+" ||
               particleName == "omega-" ||
               particleName == "sigma_c+" ||
               particleName == "sigma_c++" ||
               particleName == "sigma+" ||
               particleName == "sigma-" ||
               particleName == "tau+" ||
               particleName == "tau-" ||
               particleName == "triton" ||
               particleName == "xi_c+" ||
               particleName == "xi-"       ) 
    {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
    }
  }
    
  // Em options
  //    
  G4EmProcessOptions opt;
  opt.SetVerbose(verbose);
  opt.SetApplyCuts(true);
  
  // Physics tables
  //
  //opt.SetMinEnergy(100*eV);
  //opt.SetMaxEnergy(100*TeV);
  //opt.SetDEDXBinning(84);
  //opt.SetLambdaBinning(84);
  //opt.SetSplineFlag(true);
  opt.SetPolarAngleLimit(0.2);
  
  // Energy loss
  //
  //opt.SetLinearLossLimit(0.01);
  
  // Ionization
  //
  //opt.SetSubCutoff(true);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
