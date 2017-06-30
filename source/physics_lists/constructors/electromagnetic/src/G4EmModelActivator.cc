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
#include "G4EmParticleList.hh"

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
#include "G4eCoulombScatteringModel.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4ionIonisation.hh"
#include "G4KleinNishinaModel.hh"

#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4RayleighScattering.hh" 
#include "G4UrbanMscModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4LowEPComptonModel.hh"

#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeRayleighModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeAnnihilationModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmModelActivator::G4EmModelActivator(const G4String& emphys)
  : baseName(emphys)
{
  theParameters = G4EmParameters::Instance();

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
  const std::vector<G4String> regnamesMSC = theParameters->RegionsMsc();
  if(regnamesMSC.size() > 0)
  {
    ActivateEmOptions();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ActivateEmOptions()
{
  const std::vector<G4String> regnamesPhys = theParameters->RegionsPhysics();
  G4int nreg = regnamesPhys.size();
  if(0 == nreg) { return; }
  G4int verbose = theParameters->Verbose() - 1; 
  if(verbose > 0) {
    G4cout << "### G4EmModelActivator::ActivateEmOptions for " << nreg << " regions"
           << G4endl;
  }
  const std::vector<G4String> typesPhys = theParameters->TypesMsc();

  // start configuration of models
  const G4ParticleDefinition* elec = G4Electron::Electron();
  const G4ParticleDefinition* posi = G4Positron::Positron();
  const G4ParticleDefinition* phot = G4Gamma::Gamma();
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4EmConfigurator* em_config = man->EmConfigurator();
  G4VAtomDeexcitation* adeexc = man->AtomDeexcitation(); 
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4VEmModel* mod;

  for(G4int i=0; i<nreg; ++i) {
    G4String reg = regnamesPhys[i];
    if(verbose > 0) {
      G4cout << i << ". region <" << reg << ">; type <" << typesPhys[i] << "> " 
	     << G4endl;
    }
   
    if(baseName == typesPhys[i]) { continue; }

    if("G4EmStandard" == typesPhys[i]) {
      G4UrbanMscModel* msc = new G4UrbanMscModel();
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, 100*MeV);
      msc = new G4UrbanMscModel();
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e+", "msc", msc, reg, 0.0, 100*MeV);

    } else if("G4EmStandard_opt1" == typesPhys[i] || 
	      "G4EmStandard_opt2" == typesPhys[i]) {
      G4UrbanMscModel* msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fMinimal);
      msc->SetRangeFactor(0.2);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, 100*MeV);
      msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fMinimal);
      msc->SetRangeFactor(0.2);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e+", "msc", msc, reg, 0.0, 100*MeV);
      
    } else if("G4EmStandard_opt3" == typesPhys[i]) { 
      G4UrbanMscModel* msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e-", "msc", msc, reg);
      msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e+", "msc", msc, reg);
      msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fMinimal);
      msc->SetRangeFactor(0.2);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("proton", "msc", msc, reg);
      if(G4Threading::IsMasterThread()) {
        theParameters->SetDeexActiveRegion(reg, true, false, false);
      }
      theParameters->DefineRegParamForDeex(adeexc);
      FindOrAddProcess(phot, "Rayl");
      mod = new G4LivermoreRayleighModel();
      em_config->SetExtraEmModel("gamma", "Rayl", mod, reg);
      FindOrAddProcess(phot, "phot");
      mod = new G4LivermorePhotoElectricModel();
      em_config->SetExtraEmModel("gamma", "phot", mod, reg);
      FindOrAddProcess(phot, "compt");
      mod = new G4KleinNishinaModel();
      em_config->SetExtraEmModel("gamma", "compt", mod, reg);

    } else if("G4EmStandard_opt4" == typesPhys[i]) {
      G4UrbanMscModel* msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->SetRangeFactor(0.02);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, 100*MeV);
      msc = new G4UrbanMscModel();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->SetRangeFactor(0.02);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e+", "msc", msc, reg, 0.0, 100*MeV);
      if(G4Threading::IsMasterThread()) {
        theParameters->SetDeexActiveRegion(reg, true, false, false);
        theParameters->SetStepFunction(0.2, 10*um);
        theParameters->SetStepFunctionMuHad(0.2, 10*um);
      }
      theParameters->DefineRegParamForDeex(adeexc);
      FindOrAddProcess(phot, "Rayl");
      mod = new G4LivermoreRayleighModel();
      em_config->SetExtraEmModel("gamma", "Rayl", mod, reg);
      FindOrAddProcess(phot, "phot");
      mod = new G4LivermorePhotoElectricModel();
      em_config->SetExtraEmModel("gamma", "phot", mod, reg);
      FindOrAddProcess(phot, "compt");
      mod = new G4KleinNishinaModel();
      em_config->SetExtraEmModel("gamma", "compt", mod, reg);
      mod = new G4LowEPComptonModel();
      mod->SetHighEnergyLimit(20*MeV);
      em_config->SetExtraEmModel("gamma", "compt", mod, reg);

    } else if("G4EmStandardGS" == typesPhys[i]) {
      G4GoudsmitSaundersonMscModel* msc = new G4GoudsmitSaundersonMscModel();
      msc->SetStepLimitType(fUseSafetyPlus);
      msc->SetRangeFactor(0.1);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, 100*MeV);
      msc = new G4GoudsmitSaundersonMscModel();
      msc->SetStepLimitType(fUseSafetyPlus);
      msc->SetRangeFactor(0.1);
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e+", "msc", msc, reg, 0.0, 100*MeV);

    } else if("G4EmStandardWVI" == typesPhys[i]) {
      G4WentzelVIModel* msc = new G4WentzelVIModel();
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, 100*MeV);
      msc = new G4WentzelVIModel();
      msc->SetLocked(true);
      em_config->SetExtraEmModel("e+", "msc", msc, reg, 0.0, 100*MeV);

      FindOrAddProcess(elec, "CoulombScat");
      FindOrAddProcess(posi, "CoulombScat");
      mod = new G4eCoulombScatteringModel();
      em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, 100*MeV);
      mod = new G4eCoulombScatteringModel();
      em_config->SetExtraEmModel("e+", "CoulombScat", mod, reg, 0.0, 100*MeV);

      if(G4Threading::IsMasterThread()) {
        theParameters->SetDeexActiveRegion(regnamesPhys[i], true, false, false);
      }
      theParameters->DefineRegParamForDeex(adeexc);

    } else if("G4EmStandardSS" == typesPhys[i] && 
	      baseName != "G4EmStandard_opt3") {
      G4EmParticleList emList;
      for(const auto& particleName : emList.PartNames()) {
	G4ParticleDefinition* particle = table->FindParticle(particleName);
        if(particle && 0.0 != particle->GetPDGCharge()) {
	  FindOrAddProcess(particle, "CoulombScat");
	  G4eCoulombScatteringModel* sc = new G4eCoulombScatteringModel();
	  sc->SetPolarAngleLimit(0.0);
	  sc->SetLocked(true);
	  em_config->SetExtraEmModel(particleName, "CoulombScat", sc, reg);
          if(particleName == "mu+" || particleName == "mu-") {
	    em_config->SetExtraEmModel(particleName, "muMsc", 
				       new G4DummyModel(), reg);
	  } else {
	    em_config->SetExtraEmModel(particleName, "msc", 
				       new G4DummyModel(), reg);
	  }
	}
      }
      if(G4Threading::IsMasterThread()) {
        theParameters->SetDeexActiveRegion(regnamesPhys[i], true, true, true);
      }
      theParameters->DefineRegParamForDeex(adeexc);

    } else if("G4EmLivermore" == typesPhys[i]) {

      mod = new G4LivermorePhotoElectricModel();
      em_config->SetExtraEmModel("gamma", "phot", mod, reg);
      mod = new G4LivermoreComptonModel();
      em_config->SetExtraEmModel("gamma", "compt", mod, reg);
      mod = new G4LivermoreGammaConversionModel();
      em_config->SetExtraEmModel("gamma", "conv", mod, reg);

      FindOrAddProcess(phot, "Rayl");
      mod = new G4LivermoreRayleighModel();
      em_config->SetExtraEmModel("gamma", "Rayl", mod, reg);

      mod = new G4LivermoreIonisationModel();
      G4UniversalFluctuation* uf = new G4UniversalFluctuation();
      em_config->SetExtraEmModel("e-", "eIoni", mod, reg, 0.0, 0.1*MeV, uf);
      mod = new G4LivermoreBremsstrahlungModel();
      em_config->SetExtraEmModel("e-", "eBrem", mod, reg, 0.0, 1*GeV);

      if(G4Threading::IsMasterThread()) {
        theParameters->SetDeexActiveRegion(regnamesPhys[i], true, false, false);
      }
      theParameters->DefineRegParamForDeex(adeexc);

    } else if("G4EmPenelope" == typesPhys[i]) {
      mod = new G4PenelopePhotoElectricModel();
      em_config->SetExtraEmModel("gamma", "phot", mod, reg);
      mod = new G4PenelopeComptonModel();
      em_config->SetExtraEmModel("gamma", "compt", mod, reg);
      mod = new G4PenelopeGammaConversionModel();
      em_config->SetExtraEmModel("gamma", "conv", mod, reg);

      FindOrAddProcess(phot, "Rayl");
      mod = new G4PenelopeRayleighModel();
      em_config->SetExtraEmModel("gamma", "Rayl", mod, reg);

      mod = new G4PenelopeIonisationModel();
      G4UniversalFluctuation* uf = new G4UniversalFluctuation();
      em_config->SetExtraEmModel("e-", "eIoni", mod, reg, 0.0, 1*GeV, uf);
      mod = new G4PenelopeBremsstrahlungModel();
      em_config->SetExtraEmModel("e-", "eBrem", mod, reg, 0.0, 1*GeV);

      mod = new G4PenelopeIonisationModel();
      uf = new G4UniversalFluctuation();
      em_config->SetExtraEmModel("e+", "eIoni", mod, reg, 0.0, 1*GeV, uf);
      mod = new G4PenelopeBremsstrahlungModel();
      em_config->SetExtraEmModel("e+", "eBrem", mod, reg, 0.0, 1*GeV);
      mod = new G4PenelopeAnnihilationModel();
      em_config->SetExtraEmModel("e+", "annihil", mod, reg, 0.0, 1*GeV);

      if(G4Threading::IsMasterThread()) {
        theParameters->SetDeexActiveRegion(regnamesPhys[i], true, false, false);
      }
      theParameters->DefineRegParamForDeex(adeexc);

    } else if("G4RadioactiveDecay" == typesPhys[i]) {

      if(G4Threading::IsMasterThread()) {
	theParameters->SetFluo(true);
	theParameters->SetAuger(true);
	theParameters->SetAugerCascade(true);
	theParameters->SetDeexcitationIgnoreCut(true);
      }
      theParameters->DefineRegParamForDeex(adeexc);      

    } else {
      if(verbose > 0 && G4Threading::IsMasterThread()) {
        G4cout << "### G4EmModelActivator::ActivateEmOptions WARNING: \n"
	       << "    EM Physics configuration name <" << typesPhys[i]
	       << "> is not known - ignored" << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmModelActivator::ActivatePAI()
{
  const std::vector<G4String> regnamesPAI = theParameters->RegionsPAI();
  G4int nreg = regnamesPAI.size();
  if(0 == nreg) { return; }
  G4int verbose = theParameters->Verbose() - 1;
  if(verbose > 0) {
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

  for(G4int i = 0; i < nreg; ++i) {
    const G4ParticleDefinition* p = nullptr;
    if(particlesPAI[i] != "all") {
      p = G4ParticleTable::GetParticleTable()->FindParticle(particlesPAI[i]);
      if(!p) {
        G4cout << "### WARNING: ActivatePAI::FindParticle fails to find "
               << particlesPAI[i] << G4endl;
        continue;
      }
    }
    const G4Region* r = regionStore->GetRegion(regnamesPAI[i], false);
    if(!r) {
      G4cout << "### WARNING: ActivatePAI::GetRegion fails to find "
	     << regnamesPAI[i] << G4endl;
      continue;
    }

    G4String name = "hIoni";
    if(p == elec || p == posi)
      { name = "eIoni"; }
    else if (p == mupl || p == mumi)
      { name = "muIoni"; }
    else if (p == gion)
      { name = "ionIoni"; }

    for(itr = v.begin(); itr != v.end(); ++itr) {

      G4VEnergyLossProcess* proc = *itr;
      if(!proc->IsIonisationProcess()) { continue; }

      G4String namep = proc->GetProcessName();
      if(p) {        
	if(name != namep) { continue; }
      } else {
        if(namep != "hIoni" && namep != "muIoni" && namep != "eIoni")
	  { continue; }
      }
      G4VEmModel* em = nullptr;
      G4VEmFluctuationModel* fm = nullptr;
      if(typesPAI[i] == "PAIphoton") {
	G4PAIPhotModel* mod = new G4PAIPhotModel(p,"PAIPhotModel");
	em = mod;
	fm = mod;
      } else {
	G4PAIModel* mod = new G4PAIModel(p,"PAIModel");
	em = mod;
	fm = mod;
      }
      proc->AddEmModel(0, em, fm, r);
      if(verbose > 0) {
	G4cout << "### G4EmModelActivator: add <" << typesPAI[i]
	       << "> model for " << particlesPAI[i]
	       << " in the " << regnamesPAI[i] << G4endl;
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
                               new G4UniversalFluctuation());

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

void G4EmModelActivator::FindOrAddProcess(const G4ParticleDefinition* part, 
					  const G4String& name)
{
  G4ProcessManager* pm = part->GetProcessManager();
  G4ProcessVector* pv = pm->GetProcessList();
  G4int nproc = pm->GetProcessListLength();
  for(G4int i = 0; i<nproc; ++i) {
    if(((*pv)[i])->GetProcessName() == name) { return; }
  }
  if(name == "CoulombScat") {
    G4CoulombScattering* cs = new G4CoulombScattering();
    cs->SetEmModel(new G4DummyModel());
    pm->AddDiscreteProcess(cs);
  } else if(name == "Rayl") {
    G4RayleighScattering* rs = new G4RayleighScattering();
    rs->SetEmModel(new G4DummyModel());
    pm->AddDiscreteProcess(rs);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
