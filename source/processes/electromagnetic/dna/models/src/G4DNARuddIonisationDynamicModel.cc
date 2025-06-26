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
// Created 11.02.2025 V.Ivanchenko & M. Vologzhin
//                    on base of previous Rudd models
//
// Russian Goverment grant No 075-15-2024-667 23.08.2024
//

#include "G4DNARuddIonisationDynamicModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNAGenericIonsManager.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4NistManager.hh"

#include "G4IonTable.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACrossSectionDataSet* G4DNARuddIonisationDynamicModel::xsdata_hydrogen = nullptr;
G4DNACrossSectionDataSet* G4DNARuddIonisationDynamicModel::xsdata_helium = nullptr;
G4DNACrossSectionDataSet* G4DNARuddIonisationDynamicModel::xsdata_p = nullptr;
const std::vector<G4double>* G4DNARuddIonisationDynamicModel::fpWaterDensity = nullptr;

namespace
{
  const G4double scaleFactor = CLHEP::m*CLHEP::m;
  const G4double tolerance = 1*CLHEP::eV;
  const G4double Ry = 13.6*CLHEP::eV;

  // Following values provided by M. Dingfelder (priv. comm)
  const G4double Bj[5] = {12.60*CLHEP::eV, 14.70*CLHEP::eV, 18.40*CLHEP::eV,
                          32.20*CLHEP::eV, 539*CLHEP::eV};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationDynamicModel::G4DNARuddIonisationDynamicModel(const G4ParticleDefinition*,
                                                                 const G4String& nam)
  : G4VEmModel(nam)
{
  fGpow = G4Pow::GetInstance();
  fLowestEnergy = 100*CLHEP::eV;
  fAbsorptionEnergy = 50*CLHEP::eV;
  
  // Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  // Define default angular generator
  SetAngularDistribution(new G4DNARuddAngle());

  if (nullptr == xsdata_p) {
    isFirst = true;
    LoadData();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationDynamicModel::~G4DNARuddIonisationDynamicModel()
{  
  if (isFirst) {
    delete xsdata_p;
    delete xsdata_hydrogen;
    delete xsdata_helium;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::LoadData()
{
  // initialisation of static data once
  G4String filename = "dna/sigma_ionisation_p_rudd";
  xsdata_p = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
  xsdata_p->LoadData(filename);

  filename = "dna/sigma_ionisation_h_rudd";
  xsdata_hydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
  xsdata_hydrogen->LoadData(filename);

  filename = "dna/sigma_ionisation_he_rudd";
  xsdata_helium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
  xsdata_helium->LoadData(filename);
  
  // to avoid possible threading problem fill this vector only once
  auto water = G4NistManager::Instance()->FindMaterial("G4_WATER");
  fpWaterDensity =
    G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(water);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::Initialise(const G4ParticleDefinition* p,
						 const G4DataVector&)
{
  if (p != fParticle) { SetParticle(p); }

  // particle change object may be externally set
  if (nullptr == fParticleChangeForGamma) {
    fParticleChangeForGamma = GetParticleChangeForGamma();
  }
  const G4String& pname = p->GetParticleName();

  // the same definition of generic ion as in G4VEmProcess class
  if (p->GetParticleType() == "nucleus" && p->GetParticleSubType() == "generic") {
    if (pname != "deuteron" && pname != "triton" &&
	pname != "He3" && pname != "alpha" && pname != "alpha+" &&
	pname != "helium" && pname != "hydrogen") {
      isIon = true;
    }
  }

  // initialisation once in each thread
  if (!isInitialised) {
    isInitialised = true;
    xsdata = xsdata_p;

    if (pname == "helium") {
      isHelium = true;
      xsdata = xsdata_helium;
      slaterEffectiveCharge[0]=1.7;
      slaterEffectiveCharge[1]=1.15;
      slaterEffectiveCharge[2]=1.15;
      sCoefficient[0]=0.5;
      sCoefficient[1]=0.25;
      sCoefficient[2]=0.25;
      fLowestEnergy = 1*CLHEP::keV;
    } else if (pname == "alpha+") {
      isHelium = true;
      // The following values are provided by M. Dingfelder (priv. comm)
      slaterEffectiveCharge[0]=2.0;
      slaterEffectiveCharge[1]=2.0;
      slaterEffectiveCharge[2]=2.0;
      sCoefficient[0]=0.7;
      sCoefficient[1]=0.15;
      sCoefficient[2]=0.15;
    } else if (pname == "hydrogen") {
      xsdata = xsdata_hydrogen;
    }

    // defined stationary mode
    statCode = G4EmParameters::Instance()->DNAStationary();

    // initialise atomic de-excitation
    if (!statCode)
      fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

    // chemistry
    auto chem = G4DNAChemistryManager::Instance();
    if (chem->IsChemistryActivated()) {
      fChemistry = chem;
    }

    InitialiseIntegrator(0.1, 0.25, 1.05, 1*CLHEP::eV, 0.2*CLHEP::eV, 10*CLHEP::keV);

    if (verbose > 0) {
      G4cout << "### G4DNARuddIonisationDynamicModel::Initialise(..) "
	     << fParticle->GetParticleName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::SetParticle(const G4ParticleDefinition* p)
{
  fParticle = p;
  fMass = p->GetPDGMass();
  fMassRate = CLHEP::proton_mass_c2/fMass; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::StartTracking(G4Track* track)
{
  fTrack = track;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4DNARuddIonisationDynamicModel::CrossSectionPerVolume(const G4Material* material,
						       const G4ParticleDefinition* part,
						       G4double kinE,
						       G4double, G4double)
{
  // check if model is applicable for given material
  G4double density = (material->GetIndex() < fpWaterDensity->size())
    ? (*fpWaterDensity)[material->GetIndex()] : 0.0;
  if (0.0 == density) { return 0.0; }
  
  // check on kinetic energy (not scaled energy) to stop low-energy ion
  if (kinE < fAbsorptionEnergy) { return DBL_MAX; }

  // ion may be different
  if (fParticle != part) { SetParticle(part); }
  G4double q = fTrack->GetDynamicParticle()->GetCharge()*inveplus;

  // cross section for scaled energy
  G4double e = kinE*fMassRate;

  auto xs = xsdata;
  if (0.0 == q) { xs = isHelium ? xsdata_helium : xsdata_hydrogen; } 
  
  G4double sigma = (e > fLowestEnergy) ? xs->FindValue(e)
    : xs->FindValue(fLowestEnergy) * e / fLowestEnergy;

  sigma *= density;
  if (q > 1.5) { sigma *= q * q; }

  if (verbose > 1) {
    G4cout << "G4DNARuddIonisationDynamicModel for " << part->GetParticleName() 
           << " Ekin(keV)=" << kinE/CLHEP::keV 
           << " sigma(cm^2)=" << sigma/CLHEP::cm2 << G4endl;
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNARuddIonisationDynamicModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                   const G4MaterialCutsCouple* couple,
                                                   const G4DynamicParticle* dpart,
						   G4double, G4double)
{
  const G4ParticleDefinition* pd = dpart->GetDefinition();
  if (fParticle != pd) { SetParticle(pd); }

  // stop ion with energy below low energy limit
  G4double kinE = dpart->GetKineticEnergy();
  // ion shoud be stopped - check on kinetic energy and not scaled energy
  if (kinE <= fAbsorptionEnergy) {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopButAlive);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(kinE);
    return;
  }

  fScaledEnergy = kinE*fMassRate;
  fSelectedShell = SelectShell();
  G4double bindingEnergy = (useDNAWaterStructure)
    ? waterStructure.IonisationEnergy(fSelectedShell) : Bj[fSelectedShell];

  //Si: additional protection if tcs interpolation method is modified
  if (kinE < bindingEnergy) { return; }
  
  G4double esec = SampleElectronEnergy();
  G4double esum = 0.0;

  // sample deexcitation
  // here we assume that H2O electronic levels are the same as Oxygen.
  // this can be considered true with a rough 10% error in energy on K-shell,
  G4int Z = 8;	
  G4ThreeVector deltaDir = 
    GetAngularDistribution()->SampleDirectionForShell(dpart, esec, Z,
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
  G4double scatteredEnergy = kinE - bindingEnergy - esec;
  if(scatteredEnergy < -tolerance || exc < -tolerance) {
    G4cout << "G4DNARuddIonisationDynamicModel::SampleSecondaries: "
           << "negative final E(keV)=" << scatteredEnergy/CLHEP::keV << " Ein(keV)="
           << kinE/CLHEP::keV << "  " << pd->GetParticleName()
           << " Edelta(keV)=" << esec/CLHEP::keV << " MeV, Exc(keV)=" << exc/CLHEP::keV
	   << G4endl;
  }
  scatteredEnergy = std::max(scatteredEnergy, 0.0);
  
  // projectile
  if (!statCode) {
    fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(exc);
  } else {
    fParticleChangeForGamma->SetProposedKineticEnergy(kinE);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(kinE - scatteredEnergy);
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

G4int G4DNARuddIonisationDynamicModel::SelectShell()
{
  G4double sum = 0.0;
  G4double xs;
  for (G4int i=0; i<5; ++i) {
    auto ptr = xsdata->GetComponent(i);
    xs = (fScaledEnergy > fLowestEnergy) ? ptr->FindValue(fScaledEnergy)
      : ptr->FindValue(fLowestEnergy)*fScaledEnergy/fLowestEnergy;
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

G4double
G4DNARuddIonisationDynamicModel::MaxEnergy()
{
  // kinematic limit
  G4double tau = fScaledEnergy/CLHEP::proton_mass_c2;
  G4double gam = 1.0 + tau;
  G4double emax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0);

  // Initialisation of sampling
  G4double A1, B1, C1, D1, E1, A2, B2, C2, D2;
  if (fSelectedShell == 4) {
    //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
    A1 = 1.25;
    B1 = 0.5;
    C1 = 1.00;
    D1 = 1.00;
    E1 = 3.00;
    A2 = 1.10;
    B2 = 1.30;
    C2 = 1.00;
    D2 = 0.00;
    alphaConst = 0.66;
  } else {
    //Data For Liquid Water from Dingfelder (Protons in Water)
    A1 = 1.02;
    B1 = 82.0;
    C1 = 0.45;
    D1 = -0.80;
    E1 = 0.38;
    A2 = 1.07;
    // Value provided by M. Dingfelder (priv. comm)
    B2 = 11.6;
    C2 = 0.60;
    D2 = 0.04;
    alphaConst = 0.64;
  }
  bEnergy = Bj[fSelectedShell];
  G4double v2 = 0.25*emax/(bEnergy*gam*gam);
  v = std::sqrt(v2);
  u = Ry/bEnergy;
  wc = 4.*v2 - 2.*v - 0.25*u;

  G4double L1 = (C1 * fGpow->powA(v, D1)) / (1. + E1 * fGpow->powA(v, (D1 + 4.)));
  G4double L2 = C2 * fGpow->powA(v, D2);
  G4double H1 = (A1 * G4Log(1. + v2)) / (v2 + (B1 / v2));
  G4double H2 = (A2 / v2) + (B2 / (v2 * v2));

  F1 = L1 + H1;
  F2 = (L2 * H2) / (L2 + H2);
  return emax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNARuddIonisationDynamicModel::SampleElectronEnergy()
{
  // sampling is performed for proton projectile
  G4double emax = MaxEnergy();
  
  ComputeIntegral(0.0, emax);
  G4double e = SampleValue();
  if (verbose > 1) {
    G4cout << "G4DNARuddIonisationDynamicModel::SampleElectronEnergy: "
	   << fParticle->GetParticleName()
           << " Escaled(keV)=" << fScaledEnergy/CLHEP::keV << " Ee(keV)=" << e/CLHEP::keV
	   << G4endl;
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::ProbabilityDensityFunction(G4double e)
{
  // Shells ids are 0 1 2 3 4 (4 is k shell)
  // !!Attention, "energyTransfer" here is the energy transfered to the electron which means
  //             that the secondary kinetic energy is w = energyTransfer - bindingEnergy
  //
  //   ds            S                F1(nu) + w * F2(nu)
  //  ---- = G(k) * ----     -------------------------------------------
  //   dw            Bj       (1+w)^3 * [1 + exp{alpha * (w - wc) / nu}]
  //
  // w is the secondary electron kinetic Energy in eV
  //
  // All the other parameters can be found in Rudd's Papers
  //
  // M.Eugene Rudd, 1988, User-Friendly model for the energy distribution of
  // electrons from protons or electron collisions. Nucl. Tracks Rad. Meas.Vol 16 N0 2/3 pp 219-218
  //
  G4double w = e/bEnergy;
  G4double x = alphaConst*(w - wc)/v;
  G4double y = (x > -15.) ? 1.0 + G4Exp(x) : 1.0;

  G4double res = CorrectionFactor() * (F1 + w*F2) /
    (fGpow->powN((1. + w)/u, 3) * y);

  if (isHelium) {
    G4double energyTransfer = e + bEnergy;
    G4double Zeff = 2.0 -
      (sCoefficient[0] * S_1s(fScaledEnergy, energyTransfer, slaterEffectiveCharge[0], 1.) +
       sCoefficient[1] * S_2s(fScaledEnergy, energyTransfer, slaterEffectiveCharge[1], 2.) +
       sCoefficient[2] * S_2p(fScaledEnergy, energyTransfer, slaterEffectiveCharge[2], 2.) );

    res *= Zeff * Zeff;
  }
  return std::max(res, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::S_1s(G4double kine,
					       G4double energyTransfer,
					       G4double slaterEffCharge,
					       G4double shellNumber)
{
  // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (7)

  G4double r = Rh(kine, energyTransfer, slaterEffCharge, shellNumber);
  G4double value = 1. - G4Exp(-2 * r) * ( ( 2. * r + 2. ) * r + 1. );
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::S_2s(G4double kine,
					       G4double energyTransfer,
					       G4double slaterEffCharge,
					       G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (8)

  G4double r = Rh(kine, energyTransfer, slaterEffCharge, shellNumber);
  G4double value =
    1. - G4Exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::S_2p(G4double kine, 
					       G4double energyTransfer,
					       G4double slaterEffCharge,
					       G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (9)

  G4double r = Rh(kine, energyTransfer, slaterEffCharge, shellNumber);
  G4double value =
    1. - G4Exp(-2 * r) * (((( 2./3. * r + 4./3.) * r + 2.) * r + 2.) * r  + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::Rh(G4double ekin, G4double etrans,
					     G4double q, G4double shell)
{
  // The following values are provided by M. Dingfelder (priv. comm)
  // Dingfelder, in Chattanooga 2005 proceedings, p 4

  G4double escaled = CLHEP::electron_mass_c2/fMass * ekin;
  const G4double H = 13.60569172 * CLHEP::eV;
  G4double value = 2.0*std::sqrt(escaled / H)*q*H /(etrans*shell);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::CorrectionFactor() 
{
  // ZF Shortened
  G4double res = 1.0;
  if (fSelectedShell < 4) {
    const G4double ln10 = fGpow->logZ(10);
    G4double x = 2.0*((G4Log(fScaledEnergy/CLHEP::eV)/ln10) - 4.2);
    // The following values are provided by M. Dingfelder (priv. comm)
    res = 0.6/(1.0 + G4Exp(x)) + 0.9;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
