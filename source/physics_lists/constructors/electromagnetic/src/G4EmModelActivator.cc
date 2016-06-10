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
// $Id: G4EmModelActivator.cc 1651 2015-05-02 16:40:24Z vnivanch $
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmModelActivator
//
// Author:      V.Ivanchenko 01.06.2015
//
// Organisation:   G4AI
// Contract:       ESA contract 4000107387/12/NL/AT
// Customer:       ESA/ESTEC
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmModelActivator.hh"
#include "G4EmParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4Alpha.hh"
#include "G4ProcessManager.hh"
#include "G4DummyModel.hh"
#include "G4EmProcessSubType.hh"
#include "G4PhysicsListHelper.hh"

#include "G4MicroElecElastic.hh"
#include "G4MicroElecElasticModel.hh"

#include "G4MicroElecInelastic.hh"
#include "G4MicroElecInelasticModel.hh"

#include "G4BraggModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4LowECapture.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

// Processes and models for Geant4-DNA
#include "G4DNAGenericIonsManager.hh"

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAIonElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

#include "G4SystemOfUnits.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmModelActivator::G4EmModelActivator()
{
  theParameters = G4EmParameters::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmModelActivator::~G4EmModelActivator()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ConstructParticle()
{
  const std::vector<G4String> regnamesDNA = theParameters->RegionsDNA();
  if(regnamesDNA.size() > 0)
  {
    ConstructDNAParticles();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ConstructProcess()
{
  const std::vector<G4String> regnamesPAI = theParameters->RegionsPAI();
  if(regnamesPAI.size() > 0)
  {
    ActivatePAI();
  }
  const std::vector<G4String> regnamesME = theParameters->RegionsMicroElec();
  if(regnamesME.size() > 0)
  {
    ActivateMicroElec();
  }
  const std::vector<G4String> regnamesDNA = theParameters->RegionsDNA();
  if(regnamesDNA.size() > 0)
  {
    ActivateDNA();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ActivatePAI()
{
  const std::vector<G4String> regnamesPAI = theParameters->RegionsPAI();
  G4int nreg = regnamesPAI.size();
  if(0 == nreg) return;
  G4int verbose = theParameters->Verbose() - 1;
  if(verbose > 0)
  {
    G4cout << "### G4EmModelActivator::ActivatePAI for " << nreg << " regions"
           << G4endl;
  }
  const std::vector<G4String> particlesPAI = theParameters->ParticlesPAI();
  const std::vector<G4String> typesPAI = theParameters->TypesPAI();

  const std::vector<G4VEnergyLossProcess*>& v = G4LossTableManager::Instance()
      ->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  G4RegionStore* regionStore = G4RegionStore::GetInstance();

  const G4ParticleDefinition* elec = G4Electron::Electron();
  const G4ParticleDefinition* posi = G4Positron::Positron();
  const G4ParticleDefinition* mupl = G4MuonPlus::MuonPlus();
  const G4ParticleDefinition* mumi = G4MuonMinus::MuonMinus();
  const G4ParticleDefinition* gion = G4GenericIon::GenericIon();

  for(G4int i = 0; i < nreg; ++i)
  {
    G4bool res = true;
    const G4ParticleDefinition* p = 0;
    if(particlesPAI[i] != "all")
    {
      p = G4ParticleTable::GetParticleTable()->FindParticle(particlesPAI[i]);
      if(!p)
      {
        G4cout << "### WARNING: ActivatePAI::FindParticle fails to find "
               << particlesPAI[i] << G4endl;
        res = false;
      }
    }

    if(res)
    {
      const G4Region* r = regionStore->GetRegion(regnamesPAI[i], false);
      if(!r)
      {
        G4cout << "### WARNING: ActivatePAI::GetRegion fails to find "
               << regnamesPAI[i] << G4endl;
      }
      else
      {

        G4String name = "hIoni";
        if(p == elec || p == posi)
        { name = "eIoni";}
        else if (p == mupl || p == mumi)
        { name = "muIoni";}
        else if (p == gion)
        { name = "ionIoni";}

        for(itr = v.begin(); itr != v.end(); itr++)
        {
          G4VEnergyLossProcess* proc = *itr;
          if(proc->IsIonisationProcess())
          {
            if(p == 0 || (p != 0 && name == proc->GetProcessName()))
            {
              G4VEmModel* em = 0;
              G4VEmFluctuationModel* fm = 0;
              if(typesPAI[i] == "PAIphoton")
              {
                G4PAIPhotModel* mod = new G4PAIPhotModel(p,"PAIPhotModel");
                em = mod;
                fm = mod;
              }
              else
              {
                G4PAIModel* mod = new G4PAIModel(p,"PAIModel");
                em = mod;
                fm = mod;
              }
              proc->AddEmModel(0, em, fm, r);
              if(verbose > 0)
              {
                G4cout << "### G4EmModelActivator: add <" << typesPAI[i]
                << "> model for " << particlesPAI[i]
                << " in the " << regnamesPAI[i] << G4endl;
              }
            }
          }
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ActivateMicroElec()
{
  const std::vector<G4String> regnamesME = theParameters->RegionsMicroElec();
  G4int nreg = regnamesME.size();
  if(0 == nreg)
  {
    return;
  }
  G4int verbose = theParameters->Verbose() - 1;
  if(verbose > 0)
  {
    G4cout << "### G4EmModelActivator::ActivateMicroElec for " << nreg
           << " regions" << G4endl;
  }
  G4LossTableManager* man = G4LossTableManager::Instance();

  const G4ParticleDefinition* elec = G4Electron::Electron();
  const G4ParticleDefinition* prot = G4Proton::Proton();
  const G4ParticleDefinition* gion = G4GenericIon::GenericIon();
  G4ProcessManager* eman = elec->GetProcessManager();
  G4ProcessManager* pman = prot->GetProcessManager();
  G4ProcessManager* iman = gion->GetProcessManager();

  G4bool emsc = HasMsc(eman);

  // MicroElec elastic is not active in the world 
  G4MicroElecElastic* eElasProc =
      new G4MicroElecElastic("e-G4MicroElecElastic");
  eman->AddDiscreteProcess(eElasProc);

  G4MicroElecInelastic* eInelProc =
      new G4MicroElecInelastic("e-G4MicroElecInelastic");
  eman->AddDiscreteProcess(eInelProc);

  G4MicroElecInelastic* pInelProc =
      new G4MicroElecInelastic("p_G4MicroElecInelastic");
  pman->AddDiscreteProcess(pInelProc);

  G4MicroElecInelastic* iInelProc =
      new G4MicroElecInelastic("ion_G4MicroElecInelastic");
  iman->AddDiscreteProcess(iInelProc);

  // start configuration of models
  G4EmConfigurator* em_config = man->EmConfigurator();
  G4VEmModel* mod;

  // limits for MicroElec applicability
  G4double elowest = 16.7 * eV;
  G4double elimel = 100 * MeV;
  G4double elimin = 9 * MeV;
  G4double pmin = 50 * keV;
  G4double pmax = 99.9 * MeV;

  G4LowECapture* ecap = new G4LowECapture(elowest);
  eman->AddDiscreteProcess(ecap);

  for(G4int i = 0; i < nreg; ++i)
  {

    G4String reg = regnamesME[i];
    G4cout << "### MicroElec models are activated for G4Region " << reg
           << G4endl
           << "    Energy limits for e- elastic:    " << elowest/eV << " eV - "
           << elimel/MeV << " MeV"
           << G4endl
           << "    Energy limits for e- inelastic:  " << elowest/eV << " eV - "
           << elimin/MeV << " MeV"
           << G4endl
           << "    Energy limits for hadrons/ions:  " << pmin/MeV << " MeV - "
           << pmax/MeV << " MeV"
           << G4endl;

    // e-  
    if(emsc)
    {
      G4UrbanMscModel* msc = new G4UrbanMscModel();
      msc->SetActivationLowEnergyLimit(elimel);
      em_config->SetExtraEmModel("e-", "msc", msc, reg);
    }
    else
    {
      mod = new G4DummyModel();
      em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, elimel);
    }

    mod = new G4MicroElecElasticModel();
    em_config->SetExtraEmModel("e-",
                               "e-G4MicroElecElastic",
                               mod,
                               reg,
                               elowest,
                               elimin);

    mod = new G4MollerBhabhaModel();
    mod->SetActivationLowEnergyLimit(elimin);
    em_config->SetExtraEmModel("e-",
                               "eIoni",
                               mod,
                               reg,
                               0.0,
                               10 * TeV,
                               new G4UniversalFluctuation());

    mod = new G4MicroElecInelasticModel();
    em_config->SetExtraEmModel("e-",
                               "e-G4MicroElecInelastic",
                               mod,
                               reg,
                               elowest,
                               elimin);

    // proton
    mod = new G4BraggModel();
    mod->SetActivationHighEnergyLimit(pmin);
    em_config->SetExtraEmModel("proton",
                               "hIoni",
                               mod,
                               reg,
                               0.0,
                               2 * MeV,
                               new G4IonFluctuations());

    mod = new G4BetheBlochModel();
    mod->SetActivationLowEnergyLimit(pmax);
    em_config->SetExtraEmModel("proton",
                               "hIoni",
                               mod,
                               reg,
                               2 * MeV,
                               10 * TeV,
                               new G4UniversalFluctuation());

    mod = new G4MicroElecInelasticModel();
    em_config->SetExtraEmModel("proton",
                               "p_G4MicroElecInelastic",
                               mod,
                               reg,
                               pmin,
                               pmax);

    // ions
    mod = new G4BraggIonModel();
    mod->SetActivationHighEnergyLimit(pmin);
    em_config->SetExtraEmModel("GenericIon",
                               "ionIoni",
                               mod,
                               reg,
                               0.0,
                               2 * MeV,
                               new G4IonFluctuations());

    mod = new G4BetheBlochModel();
    mod->SetActivationLowEnergyLimit(pmax);
    em_config->SetExtraEmModel("GenericIon",
                               "ionIoni",
                               mod,
                               reg,
                               2 * MeV,
                               10 * TeV,
                               new G4IonFluctuations());

    mod = new G4MicroElecInelasticModel();
    em_config->SetExtraEmModel("GenericIon",
                               "ion_G4MicroElecInelastic",
                               mod,
                               reg,
                               pmin,
                               pmax);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ConstructDNAParticles()
{
  // this is only addition on top of any Physics Lists
  G4Alpha::Alpha();

  G4DNAGenericIonsManager * genericIonsManager =
      G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ActivateDNA()
{
  const std::vector<G4String>& regnamesDNA = theParameters->RegionsDNA();
  G4int nreg = regnamesDNA.size();
  if(0 == nreg)
  {
    return;
  }
  const std::vector<G4String> typesDNA = theParameters->TypesDNA();
  G4int verbose = theParameters->Verbose() - 1;
  if(verbose > 0)
  {
    G4cout << "### G4EmModelActivator::ActivateMicroElec for " << nreg
           << " regions" << G4endl;
  }
  G4LossTableManager* man = G4LossTableManager::Instance();

  // list of particles
  const G4ParticleDefinition* elec = G4Electron::Electron();
  const G4ParticleDefinition* prot = G4Proton::Proton();
  const G4ParticleDefinition* gion = G4GenericIon::GenericIon();

  G4DNAGenericIonsManager * genericIonsManager =
      G4DNAGenericIonsManager::Instance();
  const G4ParticleDefinition* alpha2 = G4Alpha::Alpha();
  const G4ParticleDefinition* alpha1 = genericIonsManager->GetIon("alpha+");
  const G4ParticleDefinition* alpha0 = genericIonsManager->GetIon("helium");
  const G4ParticleDefinition* h0 = genericIonsManager->GetIon("hydrogen");

  G4ProcessManager* eman = elec->GetProcessManager();
  G4ProcessManager* pman = prot->GetProcessManager();
  G4ProcessManager* iman = gion->GetProcessManager();
  G4ProcessManager* a2man = alpha2->GetProcessManager();
  G4ProcessManager* a1man = alpha1->GetProcessManager();
  G4ProcessManager* a0man = alpha0->GetProcessManager();
  G4ProcessManager* h0man = h0->GetProcessManager();

  G4bool emsc = HasMsc(eman);
  //G4bool pmsc  = HasMsc(pman);
  //G4bool a2msc = HasMsc(a2man);
  //G4bool a1msc = HasMsc(a1man);

  // alpha+ standard processes
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4ParticleDefinition* alpha11 = const_cast<G4ParticleDefinition*>(alpha1);
  ph->RegisterProcess(new G4hMultipleScattering(), alpha11);
  ph->RegisterProcess(new G4ionIonisation(), alpha11);

  // processes are defined with dummy models for the world 
  // elastic scatetring
  G4DNAElastic* theDNAeElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
  theDNAeElasticProcess->AddEmModel(0, new G4DummyModel());
  eman->AddDiscreteProcess(theDNAeElasticProcess);

  G4DNAElastic* theDNApElasticProcess = new G4DNAElastic("proton_G4DNAElastic");
  theDNApElasticProcess->AddEmModel(0, new G4DummyModel());
  pman->AddDiscreteProcess(theDNApElasticProcess);

  G4DNAElastic* theDNAa2ElasticProcess = new G4DNAElastic("alpha_G4DNAElastic");
  theDNAa2ElasticProcess->AddEmModel(0, new G4DummyModel());
  a2man->AddDiscreteProcess(theDNAa2ElasticProcess);

  G4DNAElastic* theDNAa1ElasticProcess =
      new G4DNAElastic("alpha+_G4DNAElastic");
  theDNAa1ElasticProcess->AddEmModel(0, new G4DummyModel());
  a1man->AddDiscreteProcess(theDNAa1ElasticProcess);

  G4DNAElastic* theDNAa0ElasticProcess =
      new G4DNAElastic("helium_G4DNAElastic");
  theDNAa0ElasticProcess->AddEmModel(0, new G4DummyModel());
  a0man->AddDiscreteProcess(theDNAa0ElasticProcess);

  G4DNAElastic* theDNAh0ElasticProcess =
      new G4DNAElastic("hydrogen_G4DNAElastic");
  theDNAh0ElasticProcess->AddEmModel(0, new G4DummyModel());
  h0man->AddDiscreteProcess(theDNAh0ElasticProcess);

  // excitation
  G4DNAExcitation* theDNAeExcProcess =
      new G4DNAExcitation("e-_G4DNAExcitation");
  theDNAeExcProcess->AddEmModel(0, new G4DummyModel());
  eman->AddDiscreteProcess(theDNAeExcProcess);

  G4DNAExcitation* theDNApExcProcess =
      new G4DNAExcitation("proton_G4DNAExcitation");
  theDNApExcProcess->AddEmModel(0, new G4DummyModel());
  pman->AddDiscreteProcess(theDNApExcProcess);

  G4DNAExcitation* theDNAa2ExcProcess =
      new G4DNAExcitation("alpha_G4DNAExcitation");
  theDNAa2ExcProcess->AddEmModel(0, new G4DummyModel());
  a2man->AddDiscreteProcess(theDNAa2ExcProcess);

  G4DNAExcitation* theDNAa1ExcProcess =
      new G4DNAExcitation("alpha+_G4DNAExcitation");
  theDNAa1ExcProcess->AddEmModel(0, new G4DummyModel());
  a1man->AddDiscreteProcess(theDNAa1ExcProcess);

  G4DNAExcitation* theDNAa0ExcProcess =
      new G4DNAExcitation("helium_G4DNAExcitation");
  theDNAa0ExcProcess->AddEmModel(0, new G4DummyModel());
  a0man->AddDiscreteProcess(theDNAa0ExcProcess);

  G4DNAExcitation* theDNAh0ExcProcess =
      new G4DNAExcitation("hydrogen_G4DNAExcitation");
  theDNAh0ExcProcess->AddEmModel(0, new G4DummyModel());
  h0man->AddDiscreteProcess(theDNAh0ExcProcess);

  // vibration excitation
  G4DNAVibExcitation* theDNAeVibExcProcess =
      new G4DNAVibExcitation("e-_G4DNAVibExcitation");
  theDNAeVibExcProcess->AddEmModel(0, new G4DummyModel());
  eman->AddDiscreteProcess(theDNAeVibExcProcess);

  // ionisation
  G4DNAIonisation* theDNAeIoniProcess =
      new G4DNAIonisation("e-_G4DNAIonisation");
  theDNAeIoniProcess->AddEmModel(0, new G4DummyModel());
  eman->AddDiscreteProcess(theDNAeIoniProcess);

  G4DNAIonisation* theDNApIoniProcess =
      new G4DNAIonisation("proton_G4DNAIonisation");
  theDNApIoniProcess->AddEmModel(0, new G4DummyModel());
  pman->AddDiscreteProcess(theDNApIoniProcess);

  G4DNAIonisation* theDNAa2IoniProcess =
      new G4DNAIonisation("alpha_G4DNAIonisation");
  theDNAa2IoniProcess->AddEmModel(0, new G4DummyModel());
  a2man->AddDiscreteProcess(theDNAa2IoniProcess);

  G4DNAIonisation* theDNAa1IoniProcess =
      new G4DNAIonisation("alpha+_G4DNAIonisation");
  theDNAa1IoniProcess->AddEmModel(0, new G4DummyModel());
  a1man->AddDiscreteProcess(theDNAa1IoniProcess);

  G4DNAIonisation* theDNAa0IoniProcess =
      new G4DNAIonisation("helium_G4DNAIonisation");
  theDNAa0IoniProcess->AddEmModel(0, new G4DummyModel());
  a0man->AddDiscreteProcess(theDNAa0IoniProcess);

  G4DNAIonisation* theDNAh0IoniProcess =
      new G4DNAIonisation("hydrogen_G4DNAIonisation");
  theDNAh0IoniProcess->AddEmModel(0, new G4DummyModel());
  h0man->AddDiscreteProcess(theDNAh0IoniProcess);

  G4DNAIonisation* theDNAiIoniProcess =
      new G4DNAIonisation("GenericIon_G4DNAIonisation");
  theDNAiIoniProcess->AddEmModel(0, new G4DummyModel());
  iman->AddDiscreteProcess(theDNAiIoniProcess);

  // attachment
  G4DNAAttachment* theDNAAttachProcess =
      new G4DNAAttachment("e-_G4DNAAttachment");
  theDNAAttachProcess->AddEmModel(0, new G4DummyModel());
  eman->AddDiscreteProcess(theDNAAttachProcess);

  // charge exchange
  G4DNAChargeDecrease* theDNApChargeDecreaseProcess =
      new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
  theDNApChargeDecreaseProcess->AddEmModel(0, new G4DummyModel());
  pman->AddDiscreteProcess(theDNApChargeDecreaseProcess);

  G4DNAChargeDecrease* theDNAa2ChargeDecreaseProcess =
      new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease");
  theDNAa2ChargeDecreaseProcess->AddEmModel(0, new G4DummyModel());
  a2man->AddDiscreteProcess(theDNAa2ChargeDecreaseProcess);

  G4DNAChargeDecrease* theDNAa1ChargeDecreaseProcess =
      new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease");
  theDNAa1ChargeDecreaseProcess->AddEmModel(0, new G4DummyModel());
  a1man->AddDiscreteProcess(theDNAa1ChargeDecreaseProcess);

  G4DNAChargeIncrease* theDNAa1ChargeIncreaseProcess =
      new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease");
  theDNAa1ChargeIncreaseProcess->AddEmModel(0, new G4DummyModel());
  a1man->AddDiscreteProcess(theDNAa1ChargeIncreaseProcess);

  G4DNAChargeIncrease* theDNAa0ChargeIncreaseProcess =
      new G4DNAChargeIncrease("helium_G4DNAChargeIncrease");
  theDNAa0ChargeIncreaseProcess->AddEmModel(0, new G4DummyModel());
  a0man->AddDiscreteProcess(theDNAa0ChargeIncreaseProcess);

  G4DNAChargeIncrease* theDNAh0ChargeIncreaseProcess =
      new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease");
  theDNAh0ChargeIncreaseProcess->AddEmModel(0, new G4DummyModel());
  h0man->AddDiscreteProcess(theDNAh0ChargeIncreaseProcess);

  // start configuration of models
  G4EmConfigurator* em_config = man->EmConfigurator();
  G4VEmModel* mod;

  // limits for DNA model applicability
  G4double elowest = 11 * eV;
  G4double elimel = 1 * MeV;
  G4double elimin = 1 * MeV;
  G4double elimvb = 100 * eV;
  G4double elimat = 13 * eV;
  G4double pmin = 10 * keV;
  G4double pminch = 100 * eV;
  G4double mgmin = 1 * keV;
  G4double gmmax = 500 * keV;
  G4double pmax = 100 * MeV;
  G4double ionmin = 10 * keV;
  G4double ionmax = 400 * MeV;
  G4double hmax = 100 * MeV;
  G4double gionmax = 1 * TeV;

  // low-energy capture
  G4LowECapture* ecap = new G4LowECapture(elowest);
  eman->AddDiscreteProcess(ecap);
  G4LowECapture* pcap = new G4LowECapture(pmin);
  pman->AddDiscreteProcess(pcap);
  G4LowECapture* icap = new G4LowECapture(ionmin);
  iman->AddDiscreteProcess(icap);
  G4LowECapture* a2cap = new G4LowECapture(ionmin);
  a2man->AddDiscreteProcess(a2cap);
  G4LowECapture* a1cap = new G4LowECapture(ionmin);
  a1man->AddDiscreteProcess(a1cap);
  G4LowECapture* a0cap = new G4LowECapture(ionmin);
  a0man->AddDiscreteProcess(a0cap);
  G4LowECapture* h0cap = new G4LowECapture(ionmin);
  h0man->AddDiscreteProcess(h0cap);

  // loop over regions
  for(G4int i = 0; i < nreg; ++i)
  {

    G4String reg = regnamesDNA[i];
    if(0 < verbose)
    {
      G4cout << "### DNA models are activated for G4Region " << reg << G4endl
          << "    Energy limits for e- elastic:    " << elowest/eV << " eV - "
          << elimel/MeV << " MeV" << G4endl
          << "    Energy limits for e- inelastic:  " << elowest/eV << " eV - "
          << elimin/MeV << " MeV" << G4endl
          << "    Energy limits for hadrons/ions:  " << pmin/MeV << " MeV - "
          << pmax/MeV << " MeV" << G4endl;
    }
    // e-
    if(emsc)
    {
      G4UrbanMscModel* msc = new G4UrbanMscModel();
      msc->SetActivationLowEnergyLimit(elimel);
      em_config->SetExtraEmModel("e-", "msc", msc, reg);
    }
    else
    {
      mod = new G4DummyModel();
      em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, elimel);
    }

    mod = new G4DNAScreenedRutherfordElasticModel();
    em_config->SetExtraEmModel("e-",
                               "e-_G4DNAElastic",
                               mod,
                               reg,
                               elowest,
                               elimel);

    mod = new G4MollerBhabhaModel();
    mod->SetActivationLowEnergyLimit(elimin);
    em_config->SetExtraEmModel("e-",
                               "eIoni",
                               mod,
                               reg,
                               0.0,
                               10 * TeV,
                               new G4UniversalFluctuation());

    mod = new G4DNABornIonisationModel();
    em_config->SetExtraEmModel("e-",
                               "e-_G4DNAIonisation",
                               mod,
                               reg,
                               elowest,
                               elimin);

    mod = new G4DNABornExcitationModel();
    em_config->SetExtraEmModel("e-",
                               "e-_G4DNAExcitation",
                               mod,
                               reg,
                               elowest,
                               elimin);

    mod = new G4DNASancheExcitationModel();
    em_config->SetExtraEmModel("e-",
                               "e-_G4DNAVibExcitation",
                               mod,
                               reg,
                               elowest,
                               elimvb);

    mod = new G4DNAMeltonAttachmentModel();
    em_config->SetExtraEmModel("e-",
                               "e-_G4DNAAttachment",
                               mod,
                               reg,
                               elowest,
                               elimat);

    // proton
    mod = new G4BraggModel();
    mod->SetActivationHighEnergyLimit(0.0);
    em_config->SetExtraEmModel("proton",
                               "hIoni",
                               mod,
                               reg,
                               0.0,
                               2 * MeV,
                               new G4IonFluctuations());

    mod = new G4BetheBlochModel();
    mod->SetActivationLowEnergyLimit(pmax);
    em_config->SetExtraEmModel("proton",
                               "hIoni",
                               mod,
                               reg,
                               2 * MeV,
                               10 * TeV,
                               new G4UniversalFluctuation());

    mod = new G4DNARuddIonisationModel();
    em_config->SetExtraEmModel("proton",
                               "proton_G4DNAIonisation",
                               mod,
                               reg,
                               0.0,
                               gmmax);

    mod = new G4DNABornIonisationModel();
    em_config->SetExtraEmModel("proton",
                               "proton_G4DNAIonisation",
                               mod,
                               reg,
                               gmmax,
                               pmax);

    mod = new G4DNAMillerGreenExcitationModel();
    em_config->SetExtraEmModel("proton",
                               "proton_G4DNAExcitation",
                               mod,
                               reg,
                               elowest,
                               gmmax);

    mod = new G4DNABornExcitationModel();
    em_config->SetExtraEmModel("proton",
                               "proton_G4DNAExcitation",
                               mod,
                               reg,
                               gmmax,
                               pmax);

    mod = new G4DNADingfelderChargeDecreaseModel();
    em_config->SetExtraEmModel("proton",
                               "proton_G4DNAChargeDecrease",
                               mod,
                               reg,
                               pminch,
                               pmax);

    mod = new G4DNAIonElasticModel();
    em_config->SetExtraEmModel("proton",
                               "proton_G4DNAIonElasticModel",
                               mod,
                               reg,
                               0.0,
                               1 * MeV);

    // ions
    mod = new G4BraggIonModel();
    mod->SetActivationHighEnergyLimit(0.0);
    em_config->SetExtraEmModel("GenericIon",
                               "ionIoni",
                               mod,
                               reg,
                               0.0,
                               2 * MeV,
                               new G4IonFluctuations());

    mod = new G4BetheBlochModel();
    mod->SetActivationLowEnergyLimit(gionmax);
    em_config->SetExtraEmModel("GenericIon",
                               "ionIoni",
                               mod,
                               reg,
                               2 * MeV,
                               10 * TeV,
                               new G4IonFluctuations());

    mod = new G4DNARuddIonisationExtendedModel();
    em_config->SetExtraEmModel("GenericIon",
                               "GenericIon_G4DNAIonisation",
                               mod,
                               reg,
                               0.0,
                               gionmax);

    // alpha++
    mod = new G4BraggIonModel();
    mod->SetActivationHighEnergyLimit(0.0);
    em_config->SetExtraEmModel("alpha",
                               "ionIoni",
                               mod,
                               reg,
                               0.0,
                               2 * MeV,
                               new G4IonFluctuations());

    mod = new G4BetheBlochModel();
    mod->SetActivationLowEnergyLimit(pmax);
    em_config->SetExtraEmModel("alpha",
                               "ionIoni",
                               mod,
                               reg,
                               2 * MeV,
                               10 * TeV,
                               new G4IonFluctuations());

    mod = new G4DNARuddIonisationExtendedModel();
    em_config->SetExtraEmModel("alpha",
                               "alpha_G4DNAIonisation",
                               mod,
                               reg,
                               0.0,
                               pmax);

    mod = new G4DNAMillerGreenExcitationModel();
    em_config->SetExtraEmModel("alpha",
                               "alpha_G4DNAExcitation",
                               mod,
                               reg,
                               mgmin,
                               pmax);

    mod = new G4DNADingfelderChargeDecreaseModel();
    em_config->SetExtraEmModel("alpha",
                               "alpha_G4DNAChargeDecrease",
                               mod,
                               reg,
                               mgmin,
                               ionmax);

    mod = new G4DNAIonElasticModel();
    em_config->SetExtraEmModel("alpha",
                               "alpha_G4DNAIonElasticModel",
                               mod,
                               reg,
                               0.0,
                               1 * MeV);

    // alpha+
    mod = new G4BraggIonModel();
    mod->SetActivationHighEnergyLimit(0.0);
    em_config->SetExtraEmModel("alpha+",
                               "ionIoni",
                               mod,
                               reg,
                               0.0,
                               2 * MeV,
                               new G4IonFluctuations());

    mod = new G4BetheBlochModel();
    mod->SetActivationLowEnergyLimit(pmax);
    em_config->SetExtraEmModel("alpha+",
                               "ionIoni",
                               mod,
                               reg,
                               2 * MeV,
                               10 * TeV,
                               new G4IonFluctuations());

    mod = new G4DNARuddIonisationModel();
    em_config->SetExtraEmModel("alpha+",
                               "alpha+_G4DNAIonisation",
                               mod,
                               reg,
                               0.0,
                               pmax);

    mod = new G4DNAMillerGreenExcitationModel();
    em_config->SetExtraEmModel("alpha+",
                               "alpha+_G4DNAExcitation",
                               mod,
                               reg,
                               mgmin,
                               pmax);

    mod = new G4DNADingfelderChargeDecreaseModel();
    em_config->SetExtraEmModel("alpha+",
                               "alpha+_G4DNAChargeDecrease",
                               mod,
                               reg,
                               mgmin,
                               ionmax);

    mod = new G4DNADingfelderChargeIncreaseModel();
    em_config->SetExtraEmModel("alpha+",
                               "alpha+_G4DNAChargeIncrease",
                               mod,
                               reg,
                               mgmin,
                               ionmax);

    mod = new G4DNAIonElasticModel();
    em_config->SetExtraEmModel("alpha+",
                               "alpha+_G4DNAIonElasticModel",
                               mod,
                               reg,
                               0.0,
                               1 * MeV);

    // helium
    mod = new G4DNARuddIonisationModel();
    em_config->SetExtraEmModel("helium",
                               "helium_G4DNAIonisation",
                               mod,
                               reg,
                               0.0,
                               pmax);

    mod = new G4DNAMillerGreenExcitationModel();
    em_config->SetExtraEmModel("helium",
                               "helium_G4DNAExcitation",
                               mod,
                               reg,
                               mgmin,
                               pmax);

    mod = new G4DNADingfelderChargeIncreaseModel();
    em_config->SetExtraEmModel("helium",
                               "helium_G4DNAChargeIncrease",
                               mod,
                               reg,
                               mgmin,
                               ionmax);

    mod = new G4DNAIonElasticModel();
    em_config->SetExtraEmModel("helium",
                               "helium_G4DNAIonElasticModel",
                               mod,
                               reg,
                               0.0,
                               1 * MeV);

    // hydrogen
    mod = new G4DNARuddIonisationModel();
    em_config->SetExtraEmModel("hydrogen",
                               "hydrogen_G4DNAIonisation",
                               mod,
                               reg,
                               0.0,
                               hmax);

    mod = new G4DNAMillerGreenExcitationModel();
    em_config->SetExtraEmModel("hydrogen",
                               "hydrogen_G4DNAExcitation",
                               mod,
                               reg,
                               elowest,
                               gmmax);

    mod = new G4DNADingfelderChargeIncreaseModel();
    em_config->SetExtraEmModel("hydrogen",
                               "hydrogen_G4DNAChargeIncrease",
                               mod,
                               reg,
                               pminch,
                               pmax);

    mod = new G4DNAIonElasticModel();
    em_config->SetExtraEmModel("hydrogen",
                               "hydrogen_G4DNAIonElasticModel",
                               mod,
                               reg,
                               0.0,
                               1 * MeV);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmModelActivator::HasMsc(G4ProcessManager* pm) const
{
  G4bool res = false;
  G4ProcessVector* pv = pm->GetProcessList();
  G4int nproc = pm->GetProcessListLength();
  for(G4int i = 0; i < nproc; ++i)
  {
    if(((*pv)[i])->GetProcessSubType() == fMultipleScattering)
    {
      res = true;
      break;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
