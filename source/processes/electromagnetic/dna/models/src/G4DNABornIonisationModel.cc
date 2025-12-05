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

#include "G4DNABornIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4NistManager.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNABornAngle.hh"
#include "G4DNASamplingTable.hh"
#include "G4LogLogInterpolation.hh"
#include "G4DeltaAngle.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACrossSectionDataSet* G4DNABornIonisationModel::xsdata_e = nullptr;
G4DNACrossSectionDataSet* G4DNABornIonisationModel::xsdata_p = nullptr;
G4DNASamplingTable* G4DNABornIonisationModel::sampling_e = nullptr;
G4DNASamplingTable* G4DNABornIonisationModel::sampling_p = nullptr;
const std::vector<G4double>* G4DNABornIonisationModel::fpWaterDensity = nullptr;

namespace
{
  G4double scaleFactor = (1.e-22 / 3.343) * CLHEP::m*CLHEP::m;
  G4double tolerance = 10*CLHEP::eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornIonisationModel::G4DNABornIonisationModel(const G4ParticleDefinition*,
                                                   const G4String& nam) :
  G4VEmModel(nam)
{
  SetDeexcitationFlag(true);

  // Define default angular generator
  SetAngularDistribution(new G4DNABornAngle());

  fasterCode = G4EmParameters::Instance()->DNAFast();
  
  if (nullptr == xsdata_p) {
    isFirst = true;
    LoadData();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornIonisationModel::~G4DNABornIonisationModel()
{
  if (isFirst) {
    delete xsdata_e;
    xsdata_e = nullptr;
    delete xsdata_p;
    xsdata_p = nullptr;
    delete sampling_e;
    sampling_e = nullptr;
    delete sampling_p;
    sampling_p = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel::LoadData()
{
  // initialisation of static data once
  G4String fileElectron("dna/sigma_ionisation_e_born");
  xsdata_e = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
  xsdata_e->LoadData(fileElectron);

  G4String fileProton("dna/sigma_ionisation_p_born");
  xsdata_p = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
  xsdata_p->LoadData(fileProton);

  // to avoid possible threading problem fill this vector only once
  auto water = G4NistManager::Instance()->FindMaterial("G4_WATER");
  fpWaterDensity =
    G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(water);

  G4bool verb = true;
  sampling_e = new G4DNASamplingTable(100);
  sampling_p = new G4DNASamplingTable(100);

  if (fasterCode) {
    G4String eb = "/dna/sigmadiff_cumulated_ionisation_e_born.dat";
    sampling_e->LoadData(eb, CLHEP::eV, 1.0, verb);
    G4String pb = "/dna/sigmadiff_cumulated_ionisation_p_born.dat";
    sampling_p->LoadData(pb, CLHEP::eV, 1.0, verb);
  } else {
    G4String eb = "/dna/sigmadiff_ionisation_e_born.dat";
    sampling_e->LoadData(eb, CLHEP::eV, scaleFactor, verb);
    G4String pb = "/dna/sigmadiff_ionisation_p_born.dat";
    sampling_p->LoadData(pb, CLHEP::eV, scaleFactor, verb);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel::Initialise(const G4ParticleDefinition* p,
					  const G4DataVector&)
{
  if (isInitialised) { return; }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
  
  if (p == G4Electron::Electron()) {
    fParticle = p;
    xsdata = xsdata_e;
    sampling = sampling_e;
    fLowEnergy = 8*CLHEP::eV;
    fHighEnergy = 1*CLHEP::MeV;
    feLimitEnergy = 19*CLHEP::eV;
    fAbsorptionEnergy = 6*CLHEP::eV;
    fMass = CLHEP::electron_mass_c2;
    isElectron = true;
  } else if (p == G4Proton::Proton()) {
    fParticle = p;
    xsdata = xsdata_p;
    sampling = sampling_p;
    fLowEnergy = 100*CLHEP::keV;
    fHighEnergy = 100*CLHEP::MeV;
    fpLimitEnergy = 70*CLHEP::MeV;
    fAbsorptionEnergy = 50*CLHEP::eV;
    fMass = CLHEP::proton_mass_c2;
    isElectron = false;
  } else {
    G4ExceptionDescription ed;
    ed << "Born ionisation model is used for " << p->GetParticleName();
    G4Exception("G4DNABornIonisationModel::Initialise","em0003",
		FatalException, ed, " it is not available.");
  }
  verbose = G4EmParameters::Instance()->WorkerVerbose();

  // defined stationary mode
  statCode = G4EmParameters::Instance()->DNAStationary();

  // initialise atomic de-excitation
  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

  // chemistry
  auto chem = G4DNAChemistryManager::Instance();
  if (chem->IsChemistryActivated()) {
    fChemistry = chem;
  }
  
  InitialiseIntegrator(0.1, 0.25, 1.05, 1*CLHEP::eV, 0.2*CLHEP::eV, 10*CLHEP::keV);
  if (verbose > 1) {
    G4cout << "Born ionisation model is initialized for " 
	   << fParticle->GetParticleName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel::StartTracking(G4Track* track)
{
  fTrack = track;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNABornIonisationModel::CrossSectionPerVolume(const G4Material* material,
                                                         const G4ParticleDefinition*,
							 G4double ekin, G4double, G4double)
{
  // check if model is applicable for given material
  G4double density = (material->GetIndex() < fpWaterDensity->size())
    ? (*fpWaterDensity)[material->GetIndex()] : 0.0;
  if (0.0 == density) { return 0.0; }
  
  // check on kinetic energy (not scaled energy) to stop low-energy ion
  const G4double xSecMax = 1.e+10*CLHEP::barn;
  if (ekin < fAbsorptionEnergy) { return xSecMax; }

  G4double e = std::min(ekin, fHighEnergy);  
  G4double sigma = (e > fLowEnergy) ? xsdata->FindValue(e)
    : xsdata->FindValue(fLowEnergy) * e / fLowEnergy;

  sigma *= density;
 
  // ICRU49 electronic SP scaling - ZF, SI
  if (!isElectron && spScaling && e < fpLimitEnergy) {
    const G4double A = 1.39241700556072800000e-9;
    const G4double B = -8.52610412942622630000e-2;
    sigma *= G4Exp(A*(ekin/CLHEP::eV) + B);
  }
  if (verbose > 1) {
    G4cout << "G4DNABornIonisationModel for " << fParticle->GetParticleName() 
           << " Ekin(keV)=" << ekin/CLHEP::keV 
           << " sigma(cm^2)=" << sigma/CLHEP::cm2 << G4endl;
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                  const G4MaterialCutsCouple* couple,
                                                  const G4DynamicParticle* dynParticle,
                                                  G4double, G4double)
{
  fPrimaryEnergy = dynParticle->GetKineticEnergy();
  // proton shoud be stopped - check on kinetic energy
  // electrons never have such low energy
  if (fPrimaryEnergy <= fAbsorptionEnergy) {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopButAlive);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(fPrimaryEnergy);
    return;
  }
  
  fSelectedShell = SelectShell();
  G4double bindingEnergy = waterStructure.IonisationEnergy(fSelectedShell);

  //SI: additional protection if tcs interpolation method is modified
  if (fPrimaryEnergy < bindingEnergy) { return; }

  // compute max energy
  if (isElectron) {
    fMaxEnergy = 0.5*(fPrimaryEnergy - bindingEnergy);
  } else {
    G4double tau = fPrimaryEnergy/fMass;
    fMaxEnergy = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0);
  }
  // SI: The following protection is necessary to avoid infinite loops :
  //   e- ionisation cross section has non zero partial xs at 18 eV for shell 2.
  //   e-  has zero cumulated partial xs at 18 eV for shell 2.
  //   This is due to the fact that the max allowed transfered energy is
  //   (18+10.79)/2=17.025 eV and only transfered energies strictly above this
  //   value have non zero partial cross section starting at transition energy 17.12 eV.
  if (fasterCode && isElectron && 2 == fSelectedShell && fPrimaryEnergy < feLimitEnergy) {
    do {
      fSelectedShell = SelectShell();
    } while (2 == fSelectedShell);
  }

  G4double esec = fasterCode ? SampleCumulative() : SampleDifferential();
  G4double esum = 0.0;

  // sample deexcitation
  // here we assume that H2O electronic levels are the same as Oxygen.
  // this can be considered true with a rough 10% error in energy on K-shell,
  G4int Z = 8;	
  G4ThreeVector deltaDir = 
    GetAngularDistribution()->SampleDirectionForShell(dynParticle, esec, Z,
						      fSelectedShell,
						      couple->GetMaterial());

  // SI: only atomic deexcitation from K shell is considered
  if (fAtomDeexcitation != nullptr && fSelectedShell == 4) {
    auto as = G4AtomicShellEnumerator(0);
    auto ashell = fAtomDeexcitation->GetAtomicShell(Z, as);
    fAtomDeexcitation->GenerateParticles(fvect, ashell, Z, 0, 0);

    // compute energy sum from de-excitation
    for (auto const & ptr : *fvect) {
      esum += ptr->GetKineticEnergy();
    }
  }
  // check energy balance
  // remaining excitation energy of water molecule
  G4double exc = std::max(bindingEnergy - esum, 0.0);

  // remaining projectile energy
  G4double scatteredEnergy = fPrimaryEnergy - bindingEnergy - esec;
  if (scatteredEnergy < -tolerance || exc < -tolerance) {
    G4cout << "G4DNABornIonisationModel::SampleSecondaries: "
           << "final E(keV)=" << scatteredEnergy/CLHEP::keV << " Ein(keV)="
           << fPrimaryEnergy/CLHEP::keV << "  " << fParticle->GetParticleName()
           << " Edelta(keV)=" << esec/CLHEP::keV << " MeV, Exc(keV)=" << exc/CLHEP::keV
	   << G4endl;
  }
  scatteredEnergy = std::max(scatteredEnergy, 0.0);

  // projectile
  if (!statCode) {
    fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(exc);
  } else {
    fParticleChangeForGamma->SetProposedKineticEnergy(fPrimaryEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(fPrimaryEnergy - scatteredEnergy);
  }

  // delta-electron
  auto dp = new G4DynamicParticle(G4Electron::Electron(), deltaDir, esec);
  fvect->push_back(dp);

  // create radical
  if (nullptr != fChemistry) {
    fChemistry->CreateWaterMolecule(eIonizedMolecule, fSelectedShell, fTrack);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNABornIonisationModel::SelectShell()
{
  G4double sum = 0.0;
  G4double xs;
  G4double e = std::min(fPrimaryEnergy, fHighEnergy); 
  for (G4int i=0; i<5; ++i) {
    auto ptr = xsdata->GetComponent(i);
    xs = (e > fLowEnergy) ? ptr->FindValue(e)
      : ptr->FindValue(fLowEnergy) * e/fLowEnergy;
    sum += xs;
    fTemp[i] = sum;
  }
  sum *= G4UniformRand();
  for (G4int i=0; i<5; ++i) {
    if (sum <= fTemp[i]) { return i; }
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel::SampleCumulative()
{
  G4double e = sampling->SampleCumulative(fPrimaryEnergy, fSelectedShell);
  if (verbose > 1) {
    G4cout << "G4DNABornIonisationModel::SampleCumulative: "
	   << fParticle->GetParticleName()
           << " Ekin(keV)=" << fPrimaryEnergy/CLHEP::keV
	   << " Ee(keV)=" << e/CLHEP::keV << G4endl;
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel::SampleDifferential()
{
  G4double xs = ComputeIntegral(0.0, fMaxEnergy); 
  G4double e = (xs > 0.0) ? SampleValue() : G4UniformRand()*fMaxEnergy;
  if (verbose > 1) {
    G4cout << "G4DNABornIonisationModel::SampleDifferential: "
	   << fParticle->GetParticleName()
           << " Ekin(keV)=" << fPrimaryEnergy/CLHEP::keV
	   << " Ee(keV)=" << e/CLHEP::keV << G4endl;
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel::ProbabilityDensityFunction(G4double ekin)
{
  return sampling->GetValue(fPrimaryEnergy, ekin, fSelectedShell);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
