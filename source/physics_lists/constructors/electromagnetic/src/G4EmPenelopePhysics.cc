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

#include "G4EmPenelopePhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// *** Processes and models

// gamma

#include "G4PhotoElectricEffect.hh"
#include "G4PenelopePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4PenelopeComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4PenelopeGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4PenelopeRayleighModel.hh"

#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"

// e- and e+

#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4PenelopeIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4PenelopeBremsstrahlungModel.hh"

// e+ only

#include "G4eplusAnnihilation.hh"
#include "G4PenelopeAnnihilationModel.hh"

// mu

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

// hadrons

#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ePairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4alphaIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// msc models
#include "G4UrbanMscModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
//

#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"

// particles

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

//
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmPenelopePhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmPenelopePhysics::G4EmPenelopePhysics(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmPenelope"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  param->SetMinEnergy(100*eV);
  param->SetLowestElectronEnergy(100*eV);
  param->SetNumberOfBinsPerDecade(20);
  param->SetStepFunction(0.2, 10*um);
  param->SetStepFunctionMuHad(0.2, 50*um);
  param->SetUseMottCorrection(true);  
  param->SetMscStepLimitType(fUseSafetyPlus);
  param->SetMscSkin(3);            
  param->SetMscRangeFactor(0.08);
  param->SetMuHadLateralDisplacement(true);
  param->SetFluo(true);
  param->SetMaxNIELEnergy(1*MeV);
  param->SetPIXEElectronCrossSectionModel("Penelope");
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmPenelopePhysics::~G4EmPenelopePhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmPenelopePhysics::ConstructParticle()
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

  // baryons
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

void G4EmPenelopePhysics::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4LossTableManager* man = G4LossTableManager::Instance();

  // muon & hadron bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
  G4MuPairProduction* mup = new G4MuPairProduction();
  G4hBremsstrahlung* pib = new G4hBremsstrahlung();
  G4hPairProduction* pip = new G4hPairProduction();
  G4hBremsstrahlung* kb = new G4hBremsstrahlung();
  G4hPairProduction* kp = new G4hPairProduction();
  G4hBremsstrahlung* pb = new G4hBremsstrahlung();
  G4hPairProduction* pp = new G4hPairProduction();
  G4ePairProduction* ee = new G4ePairProduction();

  // muon & hadron multiple scattering
  G4MuMultipleScattering* mumsc = new G4MuMultipleScattering();
  mumsc->SetEmModel(new G4WentzelVIModel());
  G4CoulombScattering* muss = new G4CoulombScattering();
  
  G4hMultipleScattering* pimsc = new G4hMultipleScattering();
  pimsc->SetEmModel(new G4WentzelVIModel());
  G4CoulombScattering* piss = new G4CoulombScattering();

  G4hMultipleScattering* kmsc = new G4hMultipleScattering();
  kmsc->SetEmModel(new G4WentzelVIModel());
  G4CoulombScattering* kss = new G4CoulombScattering();
  
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

  // high energy limit for e+- scattering models
  G4double highEnergyLimit = G4EmParameters::Instance()->MscEnergyLimit();

  // nuclear stopping
  G4double nielEnergyLimit = G4EmParameters::Instance()->MaxNIELEnergy();
  G4NuclearStopping* pnuc = new G4NuclearStopping();
  pnuc->SetMaxKinEnergy(nielEnergyLimit);

  //Applicability range for Penelope models
  //for higher energies, the Standard models are used   
  G4double PenelopeHighEnergyLimit = 1.0*GeV;

  // Add Penelope EM Processes
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for(const auto& particleName : partList.PartNames()) {
    G4ParticleDefinition* particle = table->FindParticle(particleName);
    if (!particle) { continue; } 
    if (particleName == "gamma") {

      //Photo-electric effect
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4PEEffectFluoModel());
      G4PenelopePhotoElectricModel* thePEPenelopeModel = new 
	G4PenelopePhotoElectricModel();   
      thePEPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, thePEPenelopeModel);

      //Compton scattering
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4PenelopeComptonModel* theComptonPenelopeModel = 
	new G4PenelopeComptonModel();
      theComptonScattering->SetEmModel(new G4KleinNishinaModel());
      theComptonPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      theComptonScattering->AddEmModel(0, theComptonPenelopeModel);

      //Gamma conversion
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4PenelopeGammaConversionModel* theGCPenelopeModel = 
	new G4PenelopeGammaConversionModel();
      theGCPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      theGammaConversion->AddEmModel(0, theGCPenelopeModel);

      //Rayleigh scattering
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4PenelopeRayleighModel* theRayleighPenelopeModel = 
	new G4PenelopeRayleighModel();
      theRayleigh->SetEmModel(theRayleighPenelopeModel);

      ph->RegisterProcess(thePhotoElectricEffect, particle);
      ph->RegisterProcess(theComptonScattering, particle);
      ph->RegisterProcess(theGammaConversion, particle);
      ph->RegisterProcess(theRayleigh, particle);

    } else if (particleName == "e-") {

      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->SetEmModel(msc1);
      msc->SetEmModel(msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm); 
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);
      
      //Ionisation
      G4eIonisation* eioni = new G4eIonisation();
      G4PenelopeIonisationModel* theIoniPenelope = 
	new G4PenelopeIonisationModel();
      theIoniPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);     
      eioni->AddEmModel(0,theIoniPenelope,new G4UniversalFluctuation());
      
      //Bremsstrahlung
      G4eBremsstrahlung* brem = new G4eBremsstrahlung();
      G4PenelopeBremsstrahlungModel* theBremPenelope = new 
	G4PenelopeBremsstrahlungModel();
      theBremPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      brem->SetEmModel(theBremPenelope);

      // register processes
      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eioni, particle);
      ph->RegisterProcess(brem, particle);
      ph->RegisterProcess(ee, particle);
      ph->RegisterProcess(ss, particle);

    } else if (particleName == "e+") {
    
      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->SetEmModel(msc1);
      msc->SetEmModel(msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm); 
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);
      
      //Ionisation
      G4eIonisation* eioni = new G4eIonisation();
      G4PenelopeIonisationModel* theIoniPenelope = 
	new G4PenelopeIonisationModel();
      theIoniPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      eioni->AddEmModel(0,theIoniPenelope,new G4UniversalFluctuation());

       //Bremsstrahlung
      G4eBremsstrahlung* brem = new G4eBremsstrahlung();
      G4PenelopeBremsstrahlungModel* theBremPenelope = new 
	G4PenelopeBremsstrahlungModel();
      theBremPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      brem->SetEmModel(theBremPenelope);
      
      //Annihilation
      G4eplusAnnihilation* anni = new G4eplusAnnihilation();
      G4PenelopeAnnihilationModel* theAnnPenelope = new 
	G4PenelopeAnnihilationModel();
      theAnnPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
      anni->AddEmModel(0, theAnnPenelope);

      // register processes
      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eioni, particle);
      ph->RegisterProcess(brem, particle);
      ph->RegisterProcess(ee, particle);
      ph->RegisterProcess(anni, particle);
      ph->RegisterProcess(ss, particle);

    } else if (particleName == "mu+" ||
               particleName == "mu-"    ) {

      G4MuIonisation* muIoni = new G4MuIonisation();

      ph->RegisterProcess(mumsc, particle);
      ph->RegisterProcess(muIoni, particle);
      ph->RegisterProcess(mub, particle);
      ph->RegisterProcess(mup, particle);
      ph->RegisterProcess(muss, particle);

    } else if (particleName == "alpha" ||
               particleName == "He3" ) {
      
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      G4ionIonisation* ionIoni = new G4ionIonisation();

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(ionIoni, particle);
      ph->RegisterProcess(pnuc, particle);

    } else if (particleName == "GenericIon") {

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);

      ph->RegisterProcess(hmsc, particle);
      ph->RegisterProcess(ionIoni, particle);
      ph->RegisterProcess(pnuc, particle);

    } else if (particleName == "pi+" ||
               particleName == "pi-" ) {

      G4hIonisation* hIoni = new G4hIonisation();

      ph->RegisterProcess(pimsc, particle);
      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(pib, particle);
      ph->RegisterProcess(pip, particle);
      ph->RegisterProcess(piss, particle);

    } else if (particleName == "kaon+" ||
               particleName == "kaon-" ) {

      G4hIonisation* hIoni = new G4hIonisation();

      ph->RegisterProcess(kmsc, particle);
      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(kb, particle);
      ph->RegisterProcess(kp, particle);
      ph->RegisterProcess(kss, particle);

    } else if (particleName == "proton" ||
	       particleName == "anti_proton") {

      G4hMultipleScattering* pmsc = new G4hMultipleScattering();
      pmsc->SetEmModel(new G4WentzelVIModel());
      G4hIonisation* hIoni = new G4hIonisation();

      ph->RegisterProcess(pmsc, particle);
      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(pb, particle);
      ph->RegisterProcess(pp, particle);
      ph->RegisterProcess(new G4CoulombScattering(), particle);

    } else if (particleName == "B+" ||
	       particleName == "B-" ||
	       particleName == "D+" ||
	       particleName == "D-" ||
	       particleName == "Ds+" ||
	       particleName == "Ds-" ||
               particleName == "anti_He3" ||
               particleName == "anti_alpha" ||
               particleName == "anti_deuteron" ||
               particleName == "anti_lambda_c+" ||
               particleName == "anti_omega-" ||
               particleName == "anti_sigma_c+" ||
               particleName == "anti_sigma_c++" ||
               particleName == "anti_sigma+" ||
               particleName == "anti_sigma-" ||
               particleName == "anti_triton" ||
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
               particleName == "xi-" ) {
      
      ph->RegisterProcess(hmsc, particle);
      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(pnuc, particle);
    }
  }
    
  // Nuclear stopping
  pnuc->SetMaxKinEnergy(MeV);
  
  // Deexcitation
  //
  G4VAtomDeexcitation* deexcitation = new G4UAtomicDeexcitation();
  man->SetAtomDeexcitation(deexcitation);

  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
