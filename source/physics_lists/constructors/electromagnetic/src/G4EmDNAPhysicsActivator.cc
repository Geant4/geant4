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
// Author V.Ivanchenko 
//

#include "G4EmDNAPhysicsActivator.hh"

#include "G4EmParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4VEnergyLossProcess.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
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
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4IonFluctuations.hh"
#include "G4LowECapture.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"
#include "G4ICRU49NuclearStoppingModel.hh"
#include "G4Generator2BS.hh"

#include "G4Threading.hh"
#include "G4EmDNABuilder.hh"
#include "G4EmUtility.hh"
#include "G4PhysListUtil.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysicsActivator::G4EmDNAPhysicsActivator(G4int ver)
    : G4VPhysicsConstructor("G4EmDNAPhysicsActivator"), verbose(ver)
{
  theParameters = G4EmParameters::Instance();
  theParameters->ActivateDNA();
  theParameters->SetFluo(true);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDNAPhysicsActivator::IsVerbose() const
{
  return (0 < verbose && G4Threading::IsMasterThread());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::ConstructParticle()
{
  G4EmDNABuilder::ConstructDNAParticles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::ConstructProcess()
{
  const std::vector<G4String>& regnamesDNA = theParameters->RegionsDNA();
  std::size_t nreg = regnamesDNA.size();
  if(0 == nreg)
  {
    return;
  }

  const std::vector<G4String>& typesDNA = theParameters->TypesDNA();
  G4bool fast = theParameters->DNAFast();
  G4bool st = theParameters->DNAStationary();

  const G4double emaxDNA = 1.*CLHEP::MeV;
  const G4double emaxIonDNA = 300.*CLHEP::MeV;
  const G4double eminBorn = 500.*CLHEP::keV;
  const G4double emax = theParameters->MaxKinEnergy();

  if(IsVerbose()) {
    G4cout << "### G4EmDNAPhysicsActivator::ConstructProcess for " << nreg
           << " regions; DNA physics type " << G4endl;
  }

  // list of particles
  G4ParticleDefinition* prot = G4Proton::Proton();
  G4ParticleDefinition* gion = G4GenericIon::GenericIon();

  G4DNAGenericIonsManager * genericIonsManager =
      G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* alpha2 = G4Alpha::Alpha();
  G4ParticleDefinition* alpha1 = genericIonsManager->GetIon("alpha+");
  G4ParticleDefinition* alpha0 = genericIonsManager->GetIon("helium");
  G4ParticleDefinition* h0 = genericIonsManager->GetIon("hydrogen");

  // loop over regions
  for(std::size_t i = 0; i < nreg; ++i)
  {
    if(IsVerbose())
    {
      G4cout << "### DNA models type " << typesDNA[i] 
	     << " are activated for G4Region " << regnamesDNA[i] << G4endl;
    }
    const G4Region* reg = G4EmUtility::FindRegion(regnamesDNA[i], verbose);
    if(nullptr == reg) { continue; }
    G4int opt = 0;
    if(typesDNA[i] == "DNA_Opt1") { 
      opt = 1;
    } else if(typesDNA[i] == "DNA_Opt2") {
      opt = 2;
    } else if(typesDNA[i] == "DNA_Opt3") {
      opt = 3;
    } else if(typesDNA[i] == "DNA_Opt4") {
      opt = 4;
    } else if(typesDNA[i] == "DNA_Opt5") {
      opt = 4;
    } else if(typesDNA[i] == "DNA_Opt6") {
      opt = 6;
    } else if(typesDNA[i] == "DNA_Opt7") {
      opt = 6;
    } else if(typesDNA[i] == "DNA_Opt8") {
      opt = 8;
    }
    DeactivateElectronProcesses(emaxDNA, emax, reg);
    G4EmDNABuilder::ConstructDNAElectronPhysics(emaxDNA, opt, fast, st, reg);
    DeactivateHadronProcesses(prot, emaxDNA, emax, reg);
    G4EmDNABuilder::ConstructDNAProtonPhysics(eminBorn, emaxIonDNA, opt, fast, st, reg);
    DeactivateIonProcesses(gion, emaxIonDNA, emax, reg);
    G4EmDNABuilder::ConstructDNAIonPhysics(emax, st, reg);
    DeactivateIonProcesses(alpha2, emaxIonDNA, emax, reg);
    G4EmDNABuilder::ConstructDNALightIonPhysics(alpha2, 2, opt, emaxIonDNA, fast, st, reg);
    DeactivateHadronProcesses(alpha1, emaxIonDNA, emax, reg);
    G4EmDNABuilder::ConstructDNALightIonPhysics(alpha1, 1, opt, emaxIonDNA, fast, st, reg);
    G4EmDNABuilder::ConstructDNALightIonPhysics(alpha0, 0, opt, emaxIonDNA, fast, st, reg);
    G4EmDNABuilder::ConstructDNALightIonPhysics(h0, 0, opt, emaxIonDNA, fast, st, reg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsActivator::DeactivateElectronProcesses(const G4double emaxDNA,
                                                          const G4double emax,
							  const G4Region* reg)
{
  if(emaxDNA >= emax) { return; }
  const G4double msclimit = 100.*CLHEP::MeV;
  G4ParticleDefinition* elec = G4Electron::Electron();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  G4VProcess* p;
  if(emaxDNA < msclimit) {
    p = G4PhysListUtil::FindProcess(elec, fMultipleScattering);
    G4VMultipleScattering* msc = dynamic_cast<G4VMultipleScattering*>(p);
    G4double elim = std::min(msclimit, emax);
    if(nullptr == msc) {
      msc = new G4eMultipleScattering();
      ph->RegisterProcess(msc, elec);
    }
    auto mod = new G4GoudsmitSaundersonMscModel();
    mod->SetActivationLowEnergyLimit(emaxDNA);
    mod->SetHighEnergyLimit(elim);
    msc->AddEmModel(-2, mod, reg);
  }

  p = G4PhysListUtil::FindProcess(elec, fIonisation);
  G4VEnergyLossProcess* ptr = dynamic_cast<G4VEnergyLossProcess*>(p);
  G4VEmFluctuationModel* fluc = nullptr;
  if(nullptr == ptr) {
    ptr = new G4eIonisation();
    ph->RegisterProcess(ptr, elec);
  }
  auto modi = new G4MollerBhabhaModel();
  modi->SetActivationLowEnergyLimit(emaxDNA);
  modi->SetHighEnergyLimit(emax);
  fluc = new G4UniversalFluctuation();
  ptr->AddEmModel(-2, modi, fluc, reg);

  p = G4PhysListUtil::FindProcess(elec, fBremsstrahlung);
  ptr = dynamic_cast<G4VEnergyLossProcess*>(p);
  if(nullptr == ptr) {
    ptr = new G4eBremsstrahlung();
    ph->RegisterProcess(ptr, elec);
  }
  auto modb = new G4SeltzerBergerModel();
  modb->SetAngularDistribution(new G4Generator2BS());
  modb->SetActivationLowEnergyLimit(emaxDNA);
  modb->SetHighEnergyLimit(emax);
  fluc = nullptr;
  ptr->AddEmModel(-2, modb, fluc, reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNAPhysicsActivator::DeactivateHadronProcesses(G4ParticleDefinition* part,
						   const G4double emaxDNA,
						   const G4double emax,
						   const G4Region* reg)
{
  if(emaxDNA >= emax) { return; }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4VProcess* p = G4PhysListUtil::FindProcess(part, fMultipleScattering);
  G4VMultipleScattering* msc = dynamic_cast<G4VMultipleScattering*>(p);
  if(nullptr == msc) {
    msc = new G4hMultipleScattering();
    ph->RegisterProcess(msc, part);  
  }
  G4VMscModel* mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(emaxDNA);
  mod->SetHighEnergyLimit(emax);
  msc->AddEmModel(-2, mod, reg);

  const G4double braggmax = 2*CLHEP::MeV;
  p = G4PhysListUtil::FindProcess(part, fIonisation);
  G4VEnergyLossProcess* ptr = dynamic_cast<G4VEnergyLossProcess*>(p);
  G4VEmFluctuationModel* fluc;
  G4VEmModel* br;
  if(part == G4GenericIon::GenericIon() || part == G4Alpha::Alpha()) {
    br = new G4BraggIonModel();
    fluc = new G4IonFluctuations();
  } else {
    br = new G4BraggModel();
    fluc = new G4UniversalFluctuation();
  }
  if(nullptr == ptr) {
    if(part == G4GenericIon::GenericIon() || part == G4Alpha::Alpha()) {
      ptr = new G4ionIonisation();
    } else {
      ptr = new G4hIonisation();
    }
    ptr->SetFluctModel(fluc);
    ph->RegisterProcess(ptr, part);  
  }
  br->SetActivationLowEnergyLimit(emaxDNA);
  br->SetHighEnergyLimit(braggmax);
  ptr->AddEmModel(-2, br, fluc, reg);
 
  auto be = new G4BetheBlochModel();
  be->SetLowEnergyLimit(braggmax);
  be->SetActivationLowEnergyLimit(braggmax);
  be->SetHighEnergyLimit(emax);
  ptr->AddEmModel(-3, be, fluc, reg);

  DeactivateNuclearStopping(part, emaxDNA, reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNAPhysicsActivator::DeactivateIonProcesses(G4ParticleDefinition* part,
						const G4double emaxDNA,
						const G4double emax,
						const G4Region* reg)
{
  if(emaxDNA >= emax) { return; }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4VProcess* p = G4PhysListUtil::FindProcess(part, fMultipleScattering);
  G4VMultipleScattering* msc = dynamic_cast<G4VMultipleScattering*>(p);
  if(nullptr == msc) {
    msc = new G4hMultipleScattering();
    ph->RegisterProcess(msc, part);  
  }
  auto mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(emaxDNA);
  mod->SetHighEnergyLimit(emax);
  msc->AddEmModel(-2, mod, reg);

  const G4double braggmax = 2*CLHEP::MeV;
  p = G4PhysListUtil::FindProcess(part, fIonisation);
  G4VEnergyLossProcess* ptr = dynamic_cast<G4VEnergyLossProcess*>(p);
  G4VEmFluctuationModel* fluc = new G4IonFluctuations();
  if(nullptr == ptr) {
    ptr = new G4ionIonisation();
    ptr->SetFluctModel(fluc);
    ph->RegisterProcess(ptr, part);  
  }
  auto br = new G4BraggIonModel();
  br->SetActivationLowEnergyLimit(emaxDNA);
  br->SetHighEnergyLimit(braggmax);
  ptr->AddEmModel(-2, br, fluc, reg);

  auto be = new G4BetheBlochModel();
  be->SetLowEnergyLimit(braggmax);
  be->SetActivationLowEnergyLimit(braggmax);
  be->SetHighEnergyLimit(emax);
  ptr->AddEmModel(-3, be, fluc, reg);

  DeactivateNuclearStopping(part, emaxDNA, reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmDNAPhysicsActivator::DeactivateNuclearStopping(const G4ParticleDefinition* part,
						   const G4double emax, 
						   const G4Region* reg)
{
  G4VProcess* p = G4PhysListUtil::FindProcess(part, fNuclearStopping);
  G4NuclearStopping* ptr = dynamic_cast<G4NuclearStopping*>(p);
  if(nullptr != ptr) {
    auto mod = new G4ICRU49NuclearStoppingModel();
    mod->SetActivationLowEnergyLimit(emax);
    ptr->AddEmModel(-2, mod, reg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
