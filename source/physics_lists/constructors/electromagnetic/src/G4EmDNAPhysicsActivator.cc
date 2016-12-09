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
// $Id: G4EmDNAPhysicsActivator.cc 97499 2016-06-03 12:02:26Z matkara $
// add elastic scattering processes of proton, hydrogen, helium, alpha+, alpha++

#include "G4EmDNAPhysicsActivator.hh"

#include "G4EmParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"

#include "G4Gamma.hh"
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

#include "G4BraggModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4LowECapture.hh"
#include "G4hMultipleScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4ionIonisation.hh"

// Processes and models for Geant4-DNA
#include "G4DNAGenericIonsManager.hh"

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
//#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNAIonElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

#include "G4SystemOfUnits.hh"
#include <vector>

#include "G4DNAChemistryManager.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"

#include "G4Threading.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysicsActivator::G4EmDNAPhysicsActivator(G4int ver)
  : G4VPhysicsConstructor("G4EmDNAPhysicsActivator"), verbose(ver)
{
  theParameters = G4EmParameters::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysicsActivator::~G4EmDNAPhysicsActivator()
{}

G4bool G4EmDNAPhysicsActivator::IsVerbose() const
{
  return (0 < verbose && G4Threading::IsMasterThread());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::ConstructParticle()
{
// bosons
  G4Gamma::Gamma();

// leptons
  G4Electron::Electron();
  G4Positron::Positron();
  
// baryons
  G4Proton::Proton();

  G4GenericIon::GenericIonDefinition();

  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
  //genericIonsManager->GetIon("carbon");
  //genericIonsManager->GetIon("nitrogen");
  //genericIonsManager->GetIon("oxygen");
  //genericIonsManager->GetIon("iron");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::ConstructProcess()
{
  const std::vector<G4String>& regnamesDNA = theParameters->RegionsDNA();
  G4int nreg = regnamesDNA.size();
  if(0 == nreg)
  {
    return;
  }
  const std::vector<G4String>& typesDNA = theParameters->TypesDNA();

  if(IsVerbose())
  {
    G4cout << "### G4EmDNAPhysicsActivator::ConstructProcess for " << nreg
           << " regions; DNA physics type " << typesDNA[0] << G4endl;
  }

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

  G4bool emsc  = HasMsc(eman);
  G4bool pmsc  = HasMsc(pman);
  G4bool a2msc = HasMsc(a2man);
  G4bool a1msc = HasMsc(a1man);
  G4bool imsc  = HasMsc(iman);

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

  // limits for DNA model applicability
  static const G4double elowest= 4. * eV;   
  static const G4double elimel = 1 * MeV;
  static const G4double pminbb = 2 * MeV;
  static const G4double pmin   = 10 * keV;
  static const G4double pmax   = 100 * MeV;
  static const G4double ionmin = 10 * keV;

  // low-energy capture
  G4LowECapture* ecap = nullptr;
  if(G4DNAChemistryManager::IsActivated() == false)
  {
    // Note: G4DNAElectronSolvation could also be used instead of G4LowECapture
    ecap = new G4LowECapture(elowest);
    //    ecap->SetVerboseLevel(1);
    eman->AddDiscreteProcess(ecap);
  }
  else
  {
    // When chemistry is activated: G4DNAElectronSolvation turns the electron
    // to a solvated electron, otherwise it kills the electron at the
    // corresponding high energy limit of the model
    G4DNAElectronSolvation* solvatation =
      new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
    solvatation->AddEmModel(0, new G4DummyModel());
    eman->AddDiscreteProcess(solvatation);
  }
  
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
    if(IsVerbose()) {
      G4cout << "### DNA models type " << typesDNA[i] 
	     << " are activated for G4Region " << reg << G4endl;
    }

    // type of DNA physics
    G4int itype = 0;

    if(0 == itype) {
      AddElectronModels0(reg, ecap, emsc, elowest, elimel);
      AddProtonModels0(reg, pmsc, elimel, pminbb, pmin, pmax);
      AddHeliumModels0(reg, a1msc, a2msc, elimel, pminbb, pmin, pmax);
      AddGenericIonModels0(reg, imsc, elimel, pminbb, pmin);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::AddElectronModels0(const G4String& reg, 
						 G4LowECapture* ecap, G4bool emsc, 
						 G4double elowest, G4double elimel)
{ 
  G4EmConfigurator* em_config = 
    G4LossTableManager::Instance()->EmConfigurator();
  G4VEmModel* mod;

  static const G4double elimin = 1 * MeV;
  static const G4double elimvb = 100 * eV;
  static const G4double elimat = 13 * eV;
  static const G4double elim1  = 10 * keV;

  G4double emax = theParameters->MaxKinEnergy();

  if(IsVerbose()) {
    G4cout << "    Energy limits for e- elastic:    " 
	   << elowest/eV << " eV - " << elimel/MeV << " MeV" << G4endl
	   << "    Energy limits for e- inelastic:  " 
	   << elowest/eV << " eV - " << elimin/MeV << " MeV" << G4endl;
  }

  if(emsc) {
    G4UrbanMscModel* msc = new G4UrbanMscModel();
    msc->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, 100*MeV);
  } else {
    mod = new G4eCoulombScatteringModel();
    mod->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
  }
    
  // cuts and solvation
  if(G4DNAChemistryManager::IsActivated()) {
    // For this condition to work, the chemistry list has to be called
    // before any standard EM physics list
    mod = new G4DNATransformElectronModel();
    mod->SetHighEnergyLimit(elowest);
    em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
			       mod, reg, 0.,elowest+1.*eV);
      
    mod = new G4DNAUeharaScreenedRutherfordElasticModel();
    mod->SetLowEnergyLimit(0.);
    mod->SetActivationLowEnergyLimit(0.);
    em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
			       mod, reg, 0.0, 7.4*eV);
  } else if(ecap) {
    ecap->AddRegion(reg);

  } else {
    mod = new G4DNAOneStepThermalizationModel();
    mod->SetHighEnergyLimit(elowest);
    em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
			       mod, reg, 0., elowest);
  }
    
  mod = new G4DNAChampionElasticModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
			     mod, reg, 7.4*eV, elimel);

  // ionisation
  mod = new G4MollerBhabhaModel();
  mod->SetActivationLowEnergyLimit(elimin);
  em_config->SetExtraEmModel("e-", "eIoni",
			     mod, reg, 0.0, emax,
			     new G4UniversalFluctuation());
    
  mod = new G4DNAEmfietzoglouIonisationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
			     mod, reg, 0.0, elim1);
    
  mod = new G4DNABornIonisationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
			     mod, reg, elim1, elimin);

  // exc
  mod = new G4DNAEmfietzoglouExcitationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
			     mod, reg, 0.0, elim1);
    
  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
			     mod, reg, elim1, elimin);

  mod = new G4DNASancheExcitationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
			     mod, reg, 0.0, elimvb);

  mod = new G4DNAMeltonAttachmentModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
			     mod, reg, 0.0, elimat);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::AddProtonModels0(const G4String& reg, 
					       G4bool pmsc, G4double elimel,
					       G4double pminbb, G4double pmin,
					       G4double pmax)
{ 
  G4EmConfigurator* em_config = 
    G4LossTableManager::Instance()->EmConfigurator();
  G4VEmModel* mod;

  static const G4double pminch = 100 * eV;
  static const G4double gmmax  = 500 * keV;
  static const G4double hmax   = 100 * MeV;

  G4double emax = theParameters->MaxKinEnergy();

  if(IsVerbose()) {
    G4cout << "    Energy limits for protons:  " 
	   << pmin/MeV << " MeV - " << pmax/MeV << " MeV" << G4endl;
  }

  // proton
  if(pmsc) {
    G4UrbanMscModel* msc = new G4UrbanMscModel();
    msc->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("proton", "msc", msc, reg, 0.0, 100*MeV);
  } else {
    mod = new G4eCoulombScatteringModel();
    mod->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("proton", "CoulombScat", mod, reg, 0.0, emax);
  }

  mod = new G4BraggModel();
  mod->SetActivationLowEnergyLimit(pminbb);
  em_config->SetExtraEmModel("proton", "hIoni",
			     mod, reg, 0.0, pminbb,
			     new G4UniversalFluctuation());

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(pmax);
  em_config->SetExtraEmModel("proton", "hIoni",
			     mod, reg, pminbb, emax,
			     new G4UniversalFluctuation());

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("proton",  "proton_G4DNAIonisation",
			     mod, reg, 0.0, gmmax);

  mod = new G4DNABornIonisationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
			     mod, reg, gmmax, pmax);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
			     mod, reg, 0.0, gmmax);

  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
			     mod, reg, gmmax, pmax);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAChargeDecrease",
			     mod, reg, pminch, pmax);

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAElasticModel",
			     mod, reg, 0.0, elimel);

  // hydrogen
  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAIonisation",
			     mod, reg, 0.0, hmax);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAExcitation",
			     mod, reg, 0.0, gmmax);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAChargeIncrease",
			     mod, reg, pminch, pmax);

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAElasticModel",
			     mod, reg, 0.0, elimel);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::AddGenericIonModels0(const G4String& reg, 
						   G4bool imsc, G4double elimel,
						   G4double pminbb, G4double pmin)
{ 
  G4EmConfigurator* em_config = 
    G4LossTableManager::Instance()->EmConfigurator();
  G4VEmModel* mod;

  static const G4double gionmax= 1 * GeV;

  G4double emax = theParameters->MaxKinEnergy();

  if(IsVerbose()) {
    G4cout << "    Energy limits for GenericIon:  " 
	   << pmin/MeV << " MeV/u - " << gionmax/MeV << " MeV/u" << G4endl;
  }

  if(imsc) {
    G4UrbanMscModel* msc = new G4UrbanMscModel();
    msc->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("proton", "ionmsc", msc, reg, 0.0, 100*MeV);
  } else {
    mod = new G4IonCoulombScatteringModel();
    mod->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("proton", "CoulombScat", mod, reg, 0.0, emax);
  }

  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(pminbb);
  em_config->SetExtraEmModel("GenericIon", "ionIoni",
			     mod, reg, 0.0, pminbb,
			     new G4IonFluctuations());

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(gionmax);
  em_config->SetExtraEmModel("GenericIon", "ionIoni",
			     mod, reg, pminbb, emax,
			     new G4IonFluctuations());

  mod = new G4DNARuddIonisationExtendedModel();
  em_config->SetExtraEmModel("GenericIon", "GenericIon_G4DNAIonisation",
			     mod, reg, 0.0, gionmax);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::AddHeliumModels0(const G4String& reg, G4bool a1msc, 
					       G4bool a2msc, G4double elimel,
					       G4double pminbb, G4double pmin,
					       G4double pmax)
{ 
  G4EmConfigurator* em_config = 
    G4LossTableManager::Instance()->EmConfigurator();
  G4VEmModel* mod;

  static const G4double mgmin  = 1 * keV;
  static const G4double hemax = 400 * MeV;

  G4double emax = theParameters->MaxKinEnergy();

  if(IsVerbose()) {
    G4cout << "    Energy limits for Helium ions:  " 
	   << pmin/MeV << " MeV - " << hemax/MeV << " MeV" << G4endl;
  }

  // alpha++
  if(a2msc) {
    G4UrbanMscModel* msc = new G4UrbanMscModel();
    msc->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("alpha++", "msc", msc, reg, 0.0, 100*MeV);
  }
  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(pminbb);
  em_config->SetExtraEmModel("alpha", "ionIoni",
			     mod, reg, 0.0, pminbb,
			     new G4IonFluctuations());
  
  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(pmax);
  em_config->SetExtraEmModel("alpha", "ionIoni",
			     mod, reg, pminbb, emax,
			     new G4IonFluctuations());

  mod = new G4DNARuddIonisationExtendedModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAIonisation",
			     mod, reg, 0.0, pmax);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAExcitation",
			     mod, reg, mgmin, pmax);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAChargeDecrease",
			     mod, reg, mgmin, hemax);

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAElasticModel",
			     mod, reg, 0.0, elimel);

  // ---
  // alpha+
  if(a1msc) {
    G4UrbanMscModel* msc = new G4UrbanMscModel();
    msc->SetActivationLowEnergyLimit(elimel);
    em_config->SetExtraEmModel("alpha+", "msc", msc, reg, 0.0, 100*MeV);
  }
  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(pminbb);
  em_config->SetExtraEmModel("alpha+", "ionIoni",
			     mod, reg, 0.0, pminbb,
			     new G4IonFluctuations());

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(pmax);
  em_config->SetExtraEmModel("alpha+", "ionIoni",
			     mod, reg, pminbb, emax,
			     new G4IonFluctuations());

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAIonisation",
			     mod, reg, 0.0, pmax);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAExcitation",
			     mod, reg, mgmin, pmax);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeDecrease",
			     mod, reg, mgmin, hemax);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeIncrease",
			     mod, reg, mgmin, hemax);

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAElasticModel",
			     mod, reg, 0.0, elimel);

  // ---
  // helium
  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAIonisation",
			     mod, reg, 0.0, pmax);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAExcitation",
			     mod, reg, mgmin, pmax);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAChargeIncrease",
			     mod, reg, mgmin, hemax);

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAElasticModel",
			     mod, reg, 0.0, elimel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDNAPhysicsActivator::HasMsc(G4ProcessManager* pm) const
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
