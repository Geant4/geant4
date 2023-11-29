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
// Modified by Z. Francis, S. Incerti to handle HZE 
// && inverse rudd function sampling 26-10-2010
//
// Rewitten by V.Ivanchenko 21.05.2023
//

#include "G4EmCorrections.hh"
#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

#include "G4IonTable.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"
#include "G4AutoLock.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACrossSectionDataSet* G4DNARuddIonisationExtendedModel::xsdata[] = {nullptr};
G4DNACrossSectionDataSet* G4DNARuddIonisationExtendedModel::xshelium = nullptr;
G4DNACrossSectionDataSet* G4DNARuddIonisationExtendedModel::xsalphaplus = nullptr;
const std::vector<G4double>* G4DNARuddIonisationExtendedModel::fpWaterDensity = nullptr;

namespace
{
  G4Mutex ionDNAMutex = G4MUTEX_INITIALIZER;
  const G4double scaleFactor = CLHEP::m*CLHEP::m;
  const G4double tolerance = 1*CLHEP::eV;
  const G4double Ry = 13.6*CLHEP::eV;
  const G4double Gj[5] = {0.99, 1.11, 1.11, 0.52, 1.};

  // Following values provided by M. Dingfelder (priv. comm)
  const G4double Bj[5] = {12.60*CLHEP::eV, 14.70*CLHEP::eV, 18.40*CLHEP::eV,
                          32.20*CLHEP::eV, 539*CLHEP::eV};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::G4DNARuddIonisationExtendedModel(const G4ParticleDefinition*,
                                                                   const G4String& nam)
  : G4VEmModel(nam)
{
  fEmCorrections = G4LossTableManager::Instance()->EmCorrections();
  fGpow = G4Pow::GetInstance();
  fLowestEnergy = 100*CLHEP::eV;
  fLimitEnergy = 1*CLHEP::keV;

  // Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  // Define default angular generator
  SetAngularDistribution(new G4DNARuddAngle());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::~G4DNARuddIonisationExtendedModel()
{  
  if(isFirst) {
    for(G4int i=0; i<RUDDZMAX; ++i) { delete xsdata[i]; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::Initialise(const G4ParticleDefinition* p,
                                                  const G4DataVector&)
{
  if(p != fParticle) { SetParticle(p); }

  // initialisation of static data once
  if(nullptr == xsdata[0]) {
    G4AutoLock l(&ionDNAMutex);
    if(nullptr == xsdata[0]) {
      isFirst = true;
      G4String filename("dna/sigma_ionisation_h_rudd");
      xsdata[0] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[0]->LoadData(filename);

      filename = "dna/sigma_ionisation_p_rudd";
      xsdata[1] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[1]->LoadData(filename);

      filename = "dna/sigma_ionisation_alphaplusplus_rudd";
      xsdata[2] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[2]->LoadData(filename);

      filename = "dna/sigma_ionisation_li_rudd";
      xsdata[3] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[3]->LoadData(filename);

      filename = "dna/sigma_ionisation_be_rudd";
      xsdata[4] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[4]->LoadData(filename);

      filename = "dna/sigma_ionisation_b_rudd";
      xsdata[5] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[5]->LoadData(filename);

      filename = "dna/sigma_ionisation_c_rudd";
      xsdata[6] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[6]->LoadData(filename);

      filename = "dna/sigma_ionisation_n_rudd";
      xsdata[7] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[7]->LoadData(filename);

      filename = "dna/sigma_ionisation_o_rudd";
      xsdata[8] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[8]->LoadData(filename);

      filename = "dna/sigma_ionisation_si_rudd";
      xsdata[14] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[14]->LoadData(filename);

      filename = "dna/sigma_ionisation_fe_rudd";
      xsdata[26] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsdata[26]->LoadData(filename);
      filename = "dna/sigma_ionisation_alphaplus_rudd";
      xsalphaplus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xsalphaplus->LoadData(filename);

      filename = "dna/sigma_ionisation_he_rudd";
      xshelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
      xshelium->LoadData(filename);
    }
    // to avoid possible threading problem fill this vector only once
    auto water = G4NistManager::Instance()->FindMaterial("G4_WATER");
    fpWaterDensity =
      G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(water);

    l.unlock();
  }

  // initialisation once in each thread
  if(nullptr == fParticleChangeForGamma) {
    fParticleChangeForGamma = GetParticleChangeForGamma();
    const G4String& pname = fParticle->GetParticleName();
    if(pname == "proton") {
      idx = 1;
      xscurrent = xsdata[1];
      fElow = fLowestEnergy;
    } else if(pname == "hydrogen") {
      idx = 0; 
      xscurrent = xsdata[0];
      fElow = fLowestEnergy;
    } else if(pname == "alpha") {
      idx = 1;
      xscurrent = xsdata[2];
      isHelium = true;
      fElow = fLimitEnergy;
    } else if(pname == "alpha+") {
      idx = 1;
      isHelium = true;
      xscurrent = xsalphaplus;
      fElow = fLimitEnergy;
      // The following values are provided by M. Dingfelder (priv. comm)
      slaterEffectiveCharge[0]=2.0;
      slaterEffectiveCharge[1]=2.0;
      slaterEffectiveCharge[2]=2.0;
      sCoefficient[0]=0.7;
      sCoefficient[1]=0.15;
      sCoefficient[2]=0.15;
    } else if(pname == "helium") {
      idx = 0; 
      isHelium = true;
      fElow = fLimitEnergy;
      xscurrent = xshelium;
      slaterEffectiveCharge[0]=1.7;
      slaterEffectiveCharge[1]=1.15;
      slaterEffectiveCharge[2]=1.15;
      sCoefficient[0]=0.5;
      sCoefficient[1]=0.25;
      sCoefficient[2]=0.25;
    } else {
      isIon = true;
    }
    // defined stationary mode
    statCode = G4EmParameters::Instance()->DNAStationary();

    // initialise atomic de-excitation
    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

    if (verbose > 0) {
      G4cout << "### G4DNARuddIonisationExtendedModel::Initialise(..) " << pname 
	     << "/n    idx=" << idx << " Amass=" << fAmass 
	     << " isIon=" << isIon << " isHelium=" << isHelium << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::SetParticle(const G4ParticleDefinition* p)
{
  fParticle = p;
  fMass = p->GetPDGMass();
  fAmass = p->GetAtomicMass();

  // for generic ions idx is dynamic, -1 means that data for the ion does not exist
  if(isIon) { 
    G4int i = p->GetAtomicNumber();
    idx = -1;
    if (i < RUDDZMAX && nullptr != xsdata[i]) {
      idx = i;
      fElow = fAmass*fLowestEnergy;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4DNARuddIonisationExtendedModel::CrossSectionPerVolume(const G4Material* material,
                                                        const G4ParticleDefinition* part,
                                                        G4double kinE,
                                                        G4double, G4double)
{
  // check if model is applicable for given material
  G4double density = (material->GetIndex() < fpWaterDensity->size())
    ? (*fpWaterDensity)[material->GetIndex()] : 0.0;
  if (0.0 == density) { return 0.0; }

  // ion may be different
  if (fParticle != part) { SetParticle(part); }

  // initilise mass rate
  fMassRate = 1.0;

  // ion shoud be stopped - check on kinetic energy and not scaled energy
  if (kinE < fLowestEnergy) { return DBL_MAX; }

  G4double sigma = 0.;
  
  // use ion table if available for given energy
  // for proton, hydrogen, alpha, alpha+, and helium no scaling to proton x-section
  if (idx == 0 || idx == 1) {
    sigma = (kinE > fElow) ? xscurrent->FindValue(kinE)
      : xscurrent->FindValue(fElow)*kinE/fElow;

    // for ions with data above limit energy
  } else if (idx > 1) {
    sigma = (kinE > fElow) ? xsdata[idx]->FindValue(kinE)
      : xsdata[idx]->FindValue(fElow)*kinE/fElow;

    // scaling from proton
  } else {
    fMassRate = CLHEP::proton_mass_c2/fMass;
    G4double e = kinE*fMassRate;
    sigma = (e > fLowestEnergy) ? xsdata[1]->FindValue(e)
      : xsdata[1]->FindValue(fLowestEnergy)*e/fLowestEnergy;
    sigma *= fEmCorrections->EffectiveChargeSquareRatio(part, material, kinE);
  }
  sigma *= density;

  if (verbose > 1) {
    G4cout << "G4DNARuddIonisationExtendedModel for " << part->GetParticleName() 
           << " Ekin(keV)=" << kinE/CLHEP::keV 
           << " sigma(cm^2)=" << sigma/CLHEP::cm2 << G4endl;
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNARuddIonisationExtendedModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                    const G4MaterialCutsCouple* couple,
                                                    const G4DynamicParticle* dpart,
                                                    G4double, G4double)
{
  const G4ParticleDefinition* pd = dpart->GetDefinition();
  if (fParticle != pd) { SetParticle(pd); }

  // stop ion with energy below low energy limit
  G4double kinE = dpart->GetKineticEnergy();
  // ion shoud be stopped - check on kinetic energy and not scaled energy
  if (kinE <= fLowestEnergy) {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopButAlive);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(kinE);
    return;
  }

  G4int shell = SelectShell(kinE);
  G4double bindingEnergy = (useDNAWaterStructure)
    ? waterStructure.IonisationEnergy(shell) : Bj[shell];

  //Si: additional protection if tcs interpolation method is modified
  if (kinE < bindingEnergy) return;

  G4double esec = SampleElectronEnergy(kinE, bindingEnergy, shell);
  G4double esum = 0.0;

  // sample deexcitation
  // here we assume that H2O electronic levels are the same as Oxygen.
  // this can be considered true with a rough 10% error in energy on K-shell,
  G4int Z = 8;	
  G4ThreeVector deltaDir = 
    GetAngularDistribution()->SampleDirectionForShell(dpart, esec, Z, shell, couple->GetMaterial());

  // SI: only atomic deexcitation from K shell is considered
  if(fAtomDeexcitation != nullptr && shell == 4) {
    auto as = G4AtomicShellEnumerator(shell);
    auto ashell = fAtomDeexcitation->GetAtomicShell(Z, as);
    fAtomDeexcitation->GenerateParticles(fvect, ashell, Z, 0, 0);

    // compute energy sum from de-excitation
    std::size_t nn = fvect->size();
    for (std::size_t i=0; i<nn; ++i) {
      esum += (*fvect)[i]->GetKineticEnergy();
    }
  }
  // check energy balance
  // remaining excitation energy of water molecule
  G4double exc = bindingEnergy - esum;

  // remaining projectile energy
  G4double scatteredEnergy = kinE - bindingEnergy - esec;
  if(scatteredEnergy < -tolerance || exc < -tolerance) {
    G4cout << "G4DNARuddIonisationExtendedModel::SampleSecondaries: "
           << "negative final E(keV)=" << scatteredEnergy/CLHEP::keV << " Ein(keV)="
           << kinE/CLHEP::keV << "  " << pd->GetParticleName()
           << " Edelta(keV)=" << esec/CLHEP::keV << " MeV, Exc(keV)=" << exc/CLHEP::keV
	   << G4endl;
  }

  // projectile
  if (!statCode) {
    fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(exc);
  } else {
    fParticleChangeForGamma->SetProposedKineticEnergy(kinE);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(kinE - scatteredEnergy);
  }

  // delta-electron
  G4DynamicParticle* dp = new G4DynamicParticle(G4Electron::Electron(), deltaDir, esec);
  fvect->push_back(dp);

  // create radical
  const G4Track* theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
  G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule, shell,
							 theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARuddIonisationExtendedModel::SelectShell(G4double e)
{
  G4double sum = 0.0;
  G4double xs;
  for(G4int i=0; i<5; ++i) {
    if (idx == 0 || idx == 1) {
      auto ptr = xscurrent->GetComponent(i);
      xs = (e > fElow) ? ptr->FindValue(e) : ptr->FindValue(fElow)*e/fElow;

    } else if (idx > 1) {
      auto ptr = xsdata[idx]->GetComponent(i);
      xs = (e > fElow) ? ptr->FindValue(e) : ptr->FindValue(fElow)*e/fElow;

    } else {
      // use scaling from proton
      auto ptr = xsdata[1]->GetComponent(i);
      G4double x = e*fMassRate;
      xs = (x >= fLowestEnergy) ? ptr->FindValue(x) 
	: ptr->FindValue(fLowestEnergy)*x/fLowestEnergy;
    }
    sum += xs;
    fTemp[i] = sum;
  }
  sum *= G4UniformRand();
  for(G4int i=0; i<5; ++i) {
    if(sum <= fTemp[i]) { return i; }
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::SampleElectronEnergy(G4double kine,
                                                                G4double eexc,
                                                                G4int shell)
{
  // kinematic limit
  G4double tau = kine/fMass;
  G4double emax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0);
  // compute cumulative probability function
  G4double step = 1*CLHEP::eV;
  G4int nn = (G4int)(emax/step);
  nn = std::max(nn, 10);
  step = emax/(G4double)nn;

  // find max probability
  G4double pmax = ProbabilityFunction(kine, 0.0, eexc, shell);
  //G4cout << "E(keV)=" << kine/keV << " emax=" << emax/keV
  //       << " pmax(0)=" << pmax << " shell=" << shell << " nn=" << nn << G4endl;

  G4double e2 = 0.0; // backup energy
  G4double e0 = 0.0; // energy with max probability
  G4double e = 0.0;
  for (G4int i=0; i<nn; ++i) {
    e += step;
    G4double prob = ProbabilityFunction(kine, e, eexc, shell);
    if (prob < pmax) {
      e2 = 2*e;
      break;
    }
    pmax = prob;
    e0 = e;
  }
  //G4cout << "         E0(keV)=" << e0/keV << " pmax=" << pmax << G4endl;
  pmax *= 1.05;
  // regression method with two regions
  G4double e1 = emax;
  G4double p1 = 0.0;
  if (2*e0 < emax) {
    e1 = e0 + 0.25*(emax - e0);
    p1 = ProbabilityFunction(kine, e1, eexc, shell);
  }
  G4double s2 = p1*(emax - e1);
  s2 /= (s2 + e1*pmax);
  G4double s1 = 1.0 - s2;

  // sampling
  G4int count = 0;
  G4double ymax, y, deltae;
  for (G4int i = 0; i<100000; ++i) {
    G4double q = G4UniformRand();
    if (q <= s1) {
      ymax = pmax;
      deltae = e1 * q / s1;
    } else {
      ymax = p1;
      deltae = e1 + (emax - e1)* (q - s1) / s2;
    }
    y = ProbabilityFunction(kine, deltae, eexc, shell);
    //G4cout << "    " << i << ".  deltae=" << deltae/CLHEP::keV 
    // << " y=" << y << " ymax=" << ymax << G4endl; 
    if (y > ymax && count < 10) {
      ++count;
      G4cout << "G4DNARuddIonisationExtendedModel::SampleElectronEnergy warning: "
	     << fParticle->GetParticleName() << " E(keV)=" << kine/CLHEP::keV
	     << " Edelta(keV)=" << deltae/CLHEP::keV 
	     << " y=" << y << " ymax=" << ymax << " n=" << i << G4endl; 
    }
    if (ymax * G4UniformRand() < y) {
      return deltae;
    }
  }
  deltae = e2;
  return deltae;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::ProbabilityFunction(G4double kine,
                                                               G4double deltae,
                                                               G4double bindingEnergy,
                                                               G4int shell)
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
  G4double A1, B1, C1, D1, E1, A2, B2, C2, D2, alphaConst;
  if (shell == 4) {
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
  G4double bEnergy = Bj[shell];
  G4double w = deltae/bEnergy;
  G4double u = Ry/bEnergy;
  G4double tau = kine/fMass;
  G4double gam = 1.0 + tau;;

  G4double v2 = 0.5*CLHEP::electron_mass_c2*tau*(tau + 2.0)/(bEnergy*gam*gam);
  G4double v = std::sqrt(v2);
  G4double wc = 4.*v2 - 2.*v - 0.25*u;

  G4double x = alphaConst*(w - wc)/v;
  G4double y = (x > -15.) ? 1.0 + G4Exp(x) : 1.0;

  G4double L1 = (C1 * fGpow->powA(v, D1)) / (1. + E1 * fGpow->powA(v, (D1 + 4.)));
  G4double L2 = C2 * fGpow->powA(v, D2);
  G4double H1 = (A1 * G4Log(1. + v2)) / (v2 + (B1 / v2));
  G4double H2 = (A2 / v2) + (B2 / (v2 * v2));

  G4double F1 = L1 + H1;
  G4double F2 = (L2 * H2) / (L2 + H2);

  G4double res = CorrectionFactor(kine, shell) * (F1 + w*F2) * Gj[shell] /
    (fGpow->powN((1. + w)/u, 3) * y);

  if(isHelium) {
    G4double energyTransfer = deltae + bindingEnergy;
    G4double Zeff = 2.0 -
      (sCoefficient[0] * S_1s(kine, energyTransfer, slaterEffectiveCharge[0], 1.) +
       sCoefficient[1] * S_2s(kine, energyTransfer, slaterEffectiveCharge[1], 2.) +
       sCoefficient[2] * S_2p(kine, energyTransfer, slaterEffectiveCharge[2], 2.) );

    res *= Zeff * Zeff;
  }
  return std::max(res, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::ComputeProbabilityFunction(
         const G4ParticleDefinition* p, G4double e, G4double deltae, G4int shell)
{
  if (fParticle != p) { SetParticle(p); }
  G4double bEnergy = (useDNAWaterStructure)
    ? waterStructure.IonisationEnergy(shell) : Bj[shell];
  return ProbabilityFunction(e, deltae, bEnergy, shell);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::S_1s(G4double kine,
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

G4double G4DNARuddIonisationExtendedModel::S_2s(G4double kine,
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

G4double G4DNARuddIonisationExtendedModel::S_2p(G4double kine, 
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

G4double G4DNARuddIonisationExtendedModel::Rh(G4double ekin, G4double etrans,
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

G4double 
G4DNARuddIonisationExtendedModel::CorrectionFactor(G4double kine, G4int shell) 
{
  // ZF Shortened
  G4double res = 1.0;
  if (shell < 4 && 0 != idx) {
    const G4double ln10 = fGpow->logZ(10);
    G4double x = 2.0*((G4Log(kine/CLHEP::eV)/ln10) - 4.2);
    // The following values are provided by M. Dingfelder (priv. comm)
    res = 0.6/(1.0 + G4Exp(x)) + 0.9;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
