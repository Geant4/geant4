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
// Geant4 class G4EmDNABuilder
//
// Author V.Ivanchenko 22.05.2020
//

#include "G4EmDNABuilder.hh"
#include "G4SystemOfUnits.hh"

// particles
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4ParticleDefinition.hh"

// utilities
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4EmBuilder.hh"
#include "G4PhysicsListHelper.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4PhysListUtil.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"

// standard processes
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

// standard models
#include "G4LivermorePhotoElectricModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4LowEWentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4Generator2BS.hh"
#include "G4BraggModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4DummyModel.hh"

// DNA models
#include "G4DNAOneStepThermalizationModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNACPA100ElasticModel.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"
#include "G4DNACPA100ExcitationModel.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"
#include "G4DNACPA100IonisationModel.hh"
#include "G4DNABornIonisationModel1.hh"
#include "G4DNAMeltonAttachmentModel.hh"
#include "G4DNAIonElasticModel.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNARuddIonisationModel.hh"
#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"
#include "G4DNARPWBAExcitationModel.hh"
#include "G4DNARPWBAIonisationModel.hh"

static const G4double lowEnergyRPWBA = 100*CLHEP::MeV;
static const G4double lowEnergyMSC = 1*CLHEP::MeV;
static const G4double lowEnergyProtonIoni = 2*CLHEP::MeV;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructDNAParticles()
{
  // standard particles
  G4EmBuilder::ConstructMinimalEmSet();

  // DNA ions
  G4DNAGenericIonsManager* genericIonsManager
    = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNABuilder::ConstructStandardEmPhysics(const G4double emin_elec,
                                           const G4double emin_proton,
                                           const G4double emin_alpha,
                                           const G4double emin_ion,
                                           const G4EmDNAMscModelType mscType,
                                           const G4bool)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();
  const G4double emax = param->MaxKinEnergy();
  G4EmBuilder::PrepareEMPhysics();

  // gamma
  G4ParticleDefinition* gamma = G4Gamma::Gamma();

  // photoelectric effect - Livermore model 
  auto thePEEffect = new G4PhotoElectricEffect();
  thePEEffect->SetEmModel(new G4LivermorePhotoElectricModel());
  ph->RegisterProcess(thePEEffect, gamma);

  // Compton scattering - Klein-Nishina
  auto theComptonScattering = new G4ComptonScattering();
  theComptonScattering->SetEmModel(new G4KleinNishinaModel());
  auto cModel = new G4LowEPComptonModel();
  cModel->SetHighEnergyLimit(20*CLHEP::MeV);
  theComptonScattering->AddEmModel(0, cModel);
  ph->RegisterProcess(theComptonScattering, gamma);

  // gamma conversion - 5D model
  auto theGammaConversion = new G4GammaConversion();
  ph->RegisterProcess(theGammaConversion, gamma);

  // Rayleigh scattering - Livermore model
  auto theRayleigh = new G4RayleighScattering();
  ph->RegisterProcess(theRayleigh, gamma);

  // electron 
  if(emin_elec < emax) {
    G4ParticleDefinition* elec = G4Electron::Electron();
    auto msc_el = new G4eMultipleScattering();
    G4VMscModel* msc_model_el;
    if(mscType == dnaWVI) {
      msc_model_el = new G4LowEWentzelVIModel();
    } else if(mscType == dnaGS) {
      msc_model_el = new G4GoudsmitSaundersonMscModel();
    } else {
      msc_model_el = new G4UrbanMscModel();
    }
    msc_model_el->SetActivationLowEnergyLimit(lowEnergyMSC);
    msc_el->SetEmModel(msc_model_el);
    ph->RegisterProcess(msc_el, elec);

    auto ioni = new G4eIonisation();
    auto mb_el = new G4MollerBhabhaModel();
    mb_el->SetActivationLowEnergyLimit(emin_elec);
    ioni->SetEmModel(mb_el);
    ph->RegisterProcess(ioni, elec);

    auto brem = new G4eBremsstrahlung();
    auto sb_el = new G4SeltzerBergerModel();
    sb_el->SetActivationLowEnergyLimit(emin_elec);
    sb_el->SetHighEnergyLimit(emax);
    sb_el->SetAngularDistribution(new G4Generator2BS());
    brem->SetEmModel(sb_el);
    ph->RegisterProcess(brem, elec);
  }

  // positron
  G4ParticleDefinition* posi = G4Positron::Positron();
  auto msc_pos = new G4eMultipleScattering();
  G4VMscModel* msc_model_pos;
  if(mscType == dnaWVI) {
    msc_model_pos = new G4LowEWentzelVIModel();
  } else if(mscType == dnaGS) {
    msc_model_pos = new G4GoudsmitSaundersonMscModel();
  } else {
    msc_model_pos = new G4UrbanMscModel();
  }
  msc_pos->SetEmModel(msc_model_pos);
  ph->RegisterProcess(msc_pos, posi);
  ph->RegisterProcess(new G4eIonisation(), posi);

  auto brem = new G4eBremsstrahlung();
  auto sb = new G4SeltzerBergerModel();
  sb->SetHighEnergyLimit(emax);
  sb->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(sb);
  ph->RegisterProcess(brem, posi);
  ph->RegisterProcess(new G4eplusAnnihilation(), posi);

  // proton
  if(emin_proton < emax) {
    G4ParticleDefinition* part = G4Proton::Proton();
    StandardHadronPhysics(part, lowEnergyMSC, emin_proton, emax,
                          mscType, false);
  }

  // GenericIon
  if(emin_ion < emax) {
    G4ParticleDefinition* ion = G4GenericIon::GenericIon();
    StandardHadronPhysics(ion, lowEnergyMSC, emin_ion, emax,
                          dnaUrban, true);
  }

  // alpha
  if(emin_alpha < emax) {
    G4ParticleDefinition* part = G4Alpha::Alpha();
    StandardHadronPhysics(part, lowEnergyMSC, emin_alpha, emax,
                          dnaUrban, true);
    
    // alpha+
    G4DNAGenericIonsManager* genericIonsManager
      = G4DNAGenericIonsManager::Instance();
    part = genericIonsManager->GetIon("alpha+");
    StandardHadronPhysics(part, lowEnergyMSC, emin_alpha, emax,
                          dnaUrban, false);
  }
  // list of main standard particles
  const std::vector<G4int> chargedParticles = {
    13, -13, 211, -211, 321, -321, -2212,
    1000010020, 1000010030, 1000020030
  };
  auto msc = new G4hMultipleScattering();
  msc->SetEmModel(new G4WentzelVIModel()); 
  G4EmBuilder::ConstructBasicEmPhysics(msc, chargedParticles);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::StandardHadronPhysics(G4ParticleDefinition* part,
					   const G4double lowELimitForMSC,
					   const G4double lowELimitForIoni,
					   const G4double maxEnergy,
					   const G4EmDNAMscModelType mscType,
					   const G4bool isIon)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4hMultipleScattering* msc = new G4hMultipleScattering();
  G4VMscModel* msc_model = nullptr;
  if(mscType == dnaWVI) {
    msc_model = new G4LowEWentzelVIModel();
  } else {
    msc_model = new G4UrbanMscModel();
  }
  msc_model->SetActivationLowEnergyLimit(lowELimitForMSC);
  msc_model->SetLowEnergyLimit(lowELimitForMSC);
  msc_model->SetHighEnergyLimit(maxEnergy);
  msc->SetEmModel(msc_model);
  ph->RegisterProcess(msc, part);

  G4VEnergyLossProcess* ioni = nullptr;
  G4VEmModel* mod1 = nullptr;
  if(isIon) {
    ioni = new G4ionIonisation();
    mod1 = new G4BraggIonModel();
  } else {
    ioni = new G4hIonisation();
    mod1 = new G4BraggModel();
  }
  G4double eth = lowEnergyProtonIoni*part->GetPDGMass()/CLHEP::proton_mass_c2;
  mod1->SetActivationLowEnergyLimit(lowELimitForIoni);
  mod1->SetHighEnergyLimit(eth);
  ioni->SetEmModel(mod1);

  G4VEmModel* mod2 = new G4BetheBlochModel();
  mod2->SetActivationLowEnergyLimit(lowELimitForIoni);
  mod2->SetLowEnergyLimit(eth);
  ioni->SetEmModel(mod2);
 
  ph->RegisterProcess(ioni, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNABuilder::ConstructDNAElectronPhysics(const G4double emaxDNA,
                                            const G4int opt,
                                            const G4bool fast,
                                            const G4bool stationary,
                                            const G4Region* reg)
{
  G4ParticleDefinition* part = G4Electron::Electron();

  // limit of the Emfietzoglou models
  G4double emaxE = 0.0;
  // limit of the elastic and solvation models
  G4double emaxT = 7.4*CLHEP::eV;
  // limit for CPA100 models
  G4double emaxCPA100 = 250*CLHEP::keV;
  if(4 == opt) {
    emaxE = 10.*CLHEP::keV;
    emaxT = 10.*CLHEP::eV;
  } else if(5 < opt) {
    emaxT = 11.*CLHEP::eV;
  }

  // *** Solvation ***
  G4DNAElectronSolvation* pSolvation = FindOrBuildElectronSolvation();
  auto therm = G4DNASolvationModelFactory::GetMacroDefinedModel();
  therm->SetHighEnergyLimit(emaxT);
  pSolvation->AddEmModel(-1, therm, reg);
 
  // *** Elastic scattering ***
  auto pElasticProcess = FindOrBuildElastic(part, "e-_G4DNAElastic");
  G4VEmModel* elast = nullptr;
  G4VEmModel* elast2 = nullptr;
  if(4 == opt) {
    elast = new G4DNAUeharaScreenedRutherfordElasticModel();
  } else if(5 < opt) {
    auto mod = new G4DNACPA100ElasticModel();
    mod->SelectStationary(stationary);
    elast = mod;
    elast2 = new G4DNAChampionElasticModel();
  } else {
    elast = new G4DNAChampionElasticModel();
  }
  elast->SetHighEnergyLimit(lowEnergyMSC);
  pElasticProcess->AddEmModel(-2, elast, reg);

  if(nullptr != elast2) {
    elast->SetHighEnergyLimit(emaxCPA100);
    elast2->SetLowEnergyLimit(emaxCPA100);
    elast2->SetHighEnergyLimit(lowEnergyMSC);
    pElasticProcess->AddEmModel(-3, elast2, reg);
  }

  // *** Excitation ***
  auto theDNAExc = FindOrBuildExcitation(part, "e-_G4DNAExcitation");
  if(emaxE > 0.0) {
    auto modE = new G4DNAEmfietzoglouExcitationModel();
    theDNAExc->AddEmModel(-1, modE, reg);
    modE->SelectStationary(stationary);
    modE->SetHighEnergyLimit(emaxE);
  }
  G4VEmModel* modB = nullptr;
  G4VEmModel* modB2 = nullptr;
  if(6 == opt) {
    auto mod = new G4DNACPA100ExcitationModel();
    mod->SelectStationary(stationary);
    modB = mod;
    auto mod1 = new G4DNABornExcitationModel();
    mod1->SelectStationary(stationary);
    modB2 = mod1;
  } else {
    auto mod = new G4DNABornExcitationModel();
    mod->SelectStationary(stationary);
    modB = mod;
  }
  modB->SetLowEnergyLimit(emaxE);
  modB->SetHighEnergyLimit(emaxDNA);
  theDNAExc->AddEmModel(-2, modB, reg);
  if(nullptr != modB2) {
    modB->SetHighEnergyLimit(emaxCPA100);
    modB2->SetLowEnergyLimit(emaxCPA100);
    modB2->SetHighEnergyLimit(emaxDNA);
    theDNAExc->AddEmModel(-3, modB2, reg);
  }

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, "e-_G4DNAIonisation");
  if(emaxE > 0.0) {
    auto modE = new G4DNAEmfietzoglouIonisationModel();
    theDNAIoni->AddEmModel(-1, modE, reg);
    modE->SelectFasterComputation(fast);
    modE->SelectStationary(stationary);
    modE->SetHighEnergyLimit(emaxE);
  }
  G4VEmModel* modI = nullptr;
  G4VEmModel* modI2 = nullptr;
  if(6 == opt) {
    auto mod = new G4DNACPA100IonisationModel();
    mod->SelectStationary(stationary);
    mod->SelectFasterComputation(fast);
    modI = mod;
    auto mod1 = new G4DNABornIonisationModel();
    mod1->SelectStationary(stationary);
    modI2 = mod1;
  } else {
    auto mod = new G4DNABornIonisationModel1();
    mod->SelectStationary(stationary);
    mod->SelectFasterComputation(fast);
    modI = mod;
  }
  modI->SetLowEnergyLimit(emaxE);
  modI->SetHighEnergyLimit(emaxDNA);
  theDNAIoni->AddEmModel(-2, modI, reg);  
  if(nullptr != modI2) {
    modI->SetHighEnergyLimit(emaxCPA100);
    modI2->SetLowEnergyLimit(emaxCPA100);
    modI2->SetHighEnergyLimit(emaxDNA);
    theDNAIoni->AddEmModel(-3, modI2, reg);
  }

  if(4 != opt && 6 != opt) {
    // *** Vibrational excitation ***
    auto theDNAVibExc = FindOrBuildVibExcitation(part, "e-_G4DNAVibExcitation");
    auto modS = new G4DNASancheExcitationModel();
    theDNAVibExc->AddEmModel(-1, modS, reg);
    modS->SelectStationary(stationary);
      
    // *** Attachment ***
    auto theDNAAttach = FindOrBuildAttachment(part, "e-_G4DNAAttachment");
    auto modM = new G4DNAMeltonAttachmentModel();
    theDNAAttach->AddEmModel(-1, modM, reg);
    modM->SelectStationary(stationary);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNABuilder::ConstructDNAProtonPhysics(const G4double e1DNA,
					  const G4double emaxIonDNA,
                                          const G4int opt,
                                          const G4bool fast,
                                          const G4bool stationary,
                                          const G4Region* reg)
{
  G4EmParameters* param = G4EmParameters::Instance();
  const G4double emax = param->MaxKinEnergy();
  G4ParticleDefinition* part = G4Proton::Proton();

  // *** Elastic scattering ***
  auto pElasticProcess = FindOrBuildElastic(part, "proton_G4DNAElastic");
  auto modE = new G4DNAIonElasticModel();
  modE->SetHighEnergyLimit(lowEnergyMSC);
  modE->SelectStationary(stationary);
  pElasticProcess->AddEmModel(-1, modE, reg);

  // *** Excitation ***
  G4double e2DNA = std::min(e1DNA, lowEnergyRPWBA);
  auto theDNAExc = FindOrBuildExcitation(part, "proton_G4DNAExcitation");
  auto modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SetHighEnergyLimit(e2DNA);
  modMGE->SelectStationary(stationary);
  theDNAExc->AddEmModel(-1, modMGE, reg);

  if(e2DNA < lowEnergyRPWBA) {
    auto modB = new G4DNABornExcitationModel();
    modB->SelectStationary(stationary);
    modB->SetLowEnergyLimit(e2DNA);
    modB->SetHighEnergyLimit(lowEnergyRPWBA);
    theDNAExc->AddEmModel(-2, modB, reg);
  }
  if(lowEnergyRPWBA < emaxIonDNA) {
    auto modC = new G4DNARPWBAExcitationModel();
    modC->SelectStationary(stationary);
    modC->SetLowEnergyLimit(lowEnergyRPWBA);
    modC->SetHighEnergyLimit(emaxIonDNA);
    theDNAExc->AddEmModel(-3, modC, reg);
  }

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, "proton_G4DNAIonisation");
  G4VEmModel* modRI = nullptr;
  if(2 == opt) {
    auto mod = new G4DNARuddIonisationExtendedModel();
    mod->SelectStationary(stationary);
    modRI = mod;
  } else {
    auto mod = new G4DNARuddIonisationModel();
    mod->SelectStationary(stationary);
    modRI = mod;
  }
  modRI->SetHighEnergyLimit(e1DNA);
  theDNAIoni->AddEmModel(-1, modRI, reg);

  if(e2DNA < lowEnergyRPWBA) {
    auto modI = new G4DNABornIonisationModel1();
    modI->SelectFasterComputation(fast);
    modI->SelectStationary(stationary);
    modI->SetLowEnergyLimit(e2DNA);
    modI->SetHighEnergyLimit(lowEnergyRPWBA);
    theDNAIoni->AddEmModel(-2, modI, reg);
  }
  if(lowEnergyRPWBA < emaxIonDNA) {
    auto modJ = new G4DNARPWBAIonisationModel();
    modJ->SelectFasterComputation(fast);
    modJ->SelectStationary(stationary);
    modJ->SetLowEnergyLimit(lowEnergyRPWBA);
    modJ->SetHighEnergyLimit(emaxIonDNA);
    theDNAIoni->AddEmModel(-3, modJ, reg);
  }

  // *** Charge decrease ***
  auto theDNAChargeDecreaseProcess = 
    FindOrBuildChargeDecrease(part, "proton_G4DNAChargeDecrease");
  auto modDCD = new G4DNADingfelderChargeDecreaseModel();
  modDCD->SelectStationary(stationary);
  modDCD->SetLowEnergyLimit(0.0);
  modDCD->SetHighEnergyLimit(emax);
  theDNAChargeDecreaseProcess->AddEmModel(-1, modDCD, reg);

  FindOrBuildCapture(0.1*CLHEP::keV, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNABuilder::ConstructDNAIonPhysics(const G4double emaxIonDNA,
                                       const G4bool stationary,
                                       const G4Region* reg)
{
  G4ParticleDefinition* part = G4GenericIon::GenericIon();

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, "GenericIon_G4DNAIonisation");
  auto mod = new G4DNARuddIonisationExtendedModel();
  mod->SelectStationary(stationary);
  mod->SetHighEnergyLimit(emaxIonDNA);
  theDNAIoni->AddEmModel(-1, mod, reg);

  FindOrBuildCapture(25.0*CLHEP::keV, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4EmDNABuilder::ConstructDNALightIonPhysics(G4ParticleDefinition* part,
                                            const G4int charge,
                                            const G4int opt,
                                            const G4double emaxIonDNA,
					    const G4bool,
                                            const G4bool stationary,
                                            const G4Region* reg)
{
  G4EmParameters* param = G4EmParameters::Instance();
  const G4double emax = param->MaxKinEnergy();
  const G4String& name = part->GetParticleName();

  // *** Elastic ***
  auto theDNAElastic = FindOrBuildElastic(part, name + "_G4DNAElastic");
  auto modEI = new G4DNAIonElasticModel();
  modEI->SelectStationary(stationary);
  modEI->SetHighEnergyLimit(lowEnergyMSC);
  theDNAElastic->AddEmModel(-1, modEI, reg);

  // *** Excitation ***
  auto theDNAExc = FindOrBuildExcitation(part, name + "_G4DNAExcitation");
  auto modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SelectStationary(stationary);
  modMGE->SetLowEnergyLimit(0.0);
  modMGE->SetHighEnergyLimit(emaxIonDNA);
  theDNAExc->AddEmModel(-1, modMGE, reg);

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, name + "_G4DNAIonisation");
  G4VEmModel* modRI = nullptr;
  if(2 == opt) {
    auto mod = new G4DNARuddIonisationExtendedModel();
    mod->SelectStationary(stationary);
    modRI = mod;
  } else {
    auto mod = new G4DNARuddIonisationModel();
    mod->SelectStationary(stationary);
    modRI = mod;
  }
  modRI->SetHighEnergyLimit(emaxIonDNA);
  theDNAIoni->AddEmModel(-1, modRI, reg);

  // *** Charge increase ***
  if(2 > charge) {
    auto theDNAChargeIncrease = 
      FindOrBuildChargeIncrease(part, name + "_G4DNAChargeIncrease");
    auto modDCI = new G4DNADingfelderChargeIncreaseModel();
    modDCI->SelectStationary(stationary);
    modDCI->SetLowEnergyLimit(0.0);
    modDCI->SetHighEnergyLimit(emax);
    theDNAChargeIncrease->AddEmModel(-1, modDCI, reg);
  }

  // *** Charge decrease ***
  if(0 < charge) {
    auto theDNAChargeDecrease = 
      FindOrBuildChargeDecrease(part, name + "_G4DNAChargeDecrease");
    auto modDCD = new G4DNADingfelderChargeDecreaseModel();
    modDCD->SelectStationary(stationary);
    modDCD->SetLowEnergyLimit(0.0);
    modDCD->SetHighEnergyLimit(emax);
    theDNAChargeDecrease->AddEmModel(-1, modDCD, reg);
  }
  FindOrBuildCapture(1*CLHEP::keV, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAElectronSolvation* G4EmDNABuilder::FindOrBuildElectronSolvation()
{
  auto elec = G4Electron::Electron();
  auto* p = G4PhysListUtil::FindProcess(elec, fLowEnergyElectronSolvation);
  G4DNAElectronSolvation* ptr = dynamic_cast<G4DNAElectronSolvation*>(p);
  if(nullptr == ptr) {
    ptr = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, elec);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAElastic* 
G4EmDNABuilder::FindOrBuildElastic(G4ParticleDefinition* part,
                                   const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyElastic);
  G4DNAElastic* ptr = dynamic_cast<G4DNAElastic*>(p);
  if(nullptr == ptr) {
    ptr = new G4DNAElastic(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAExcitation* 
G4EmDNABuilder::FindOrBuildExcitation(G4ParticleDefinition* part,
                                      const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyExcitation);
  G4DNAExcitation* ptr = dynamic_cast<G4DNAExcitation*>(p);
  if(nullptr == ptr) { 
    ptr = new G4DNAExcitation(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAVibExcitation* 
G4EmDNABuilder::FindOrBuildVibExcitation(G4ParticleDefinition* part,
                                         const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyVibrationalExcitation);
  G4DNAVibExcitation* ptr = dynamic_cast<G4DNAVibExcitation*>(p);
  if(nullptr == ptr) { 
    ptr = new G4DNAVibExcitation(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAIonisation* 
G4EmDNABuilder::FindOrBuildIonisation(G4ParticleDefinition* part,
                                      const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyIonisation);
  G4DNAIonisation* ptr = dynamic_cast<G4DNAIonisation*>(p);
  if(nullptr == ptr) { 
    ptr = new G4DNAIonisation(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAAttachment* 
G4EmDNABuilder::FindOrBuildAttachment(G4ParticleDefinition* part,
                                      const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyAttachment);
  G4DNAAttachment* ptr = dynamic_cast<G4DNAAttachment*>(p);
  if(nullptr == ptr) { 
    ptr = new G4DNAAttachment(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAChargeDecrease*
G4EmDNABuilder::FindOrBuildChargeDecrease(G4ParticleDefinition* part,
                                          const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyChargeDecrease);
  G4DNAChargeDecrease* ptr = dynamic_cast<G4DNAChargeDecrease*>(p);
  if(nullptr == ptr) { 
    ptr = new G4DNAChargeDecrease(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAChargeIncrease*
G4EmDNABuilder::FindOrBuildChargeIncrease(G4ParticleDefinition* part,
                                          const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyChargeIncrease);
  G4DNAChargeIncrease* ptr = dynamic_cast<G4DNAChargeIncrease*>(p);
  if(nullptr == ptr) { 
    ptr = new G4DNAChargeIncrease(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LowECapture*
G4EmDNABuilder::FindOrBuildCapture(const G4double elim, G4ParticleDefinition* part)
{
  auto p = G4PhysListUtil::FindProcess(part, 66);
  G4LowECapture* ptr = dynamic_cast<G4LowECapture*>(p);
  if(nullptr == ptr) { 
    ptr = new G4LowECapture(elim);
    auto mng = part->GetProcessManager();
    mng->AddDiscreteProcess(ptr);
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
