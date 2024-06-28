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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACrossSectionDataSet* G4DNARuddIonisationExtendedModel::xsdata[] = {nullptr};
G4DNACrossSectionDataSet* G4DNARuddIonisationExtendedModel::xshelium = nullptr;
G4DNACrossSectionDataSet* G4DNARuddIonisationExtendedModel::xsalphaplus = nullptr;
const std::vector<G4double>* G4DNARuddIonisationExtendedModel::fpWaterDensity = nullptr;

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

  if (nullptr == xshelium) { LoadData(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::~G4DNARuddIonisationExtendedModel()
{  
  if(isFirst) {
    for(auto & i : xsdata) { delete i; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::LoadData()
{
  // initialisation of static data once
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

  // to avoid possible threading problem fill this vector only once
  auto water = G4NistManager::Instance()->FindMaterial("G4_WATER");
  fpWaterDensity =
    G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(water);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::Initialise(const G4ParticleDefinition* p,
                                                  const G4DataVector&)
{
  if (p != fParticle) { SetParticle(p); }

  // particle change object may be externally set
  if (nullptr == fParticleChangeForGamma) {
    fParticleChangeForGamma = GetParticleChangeForGamma();
  }

  // initialisation once in each thread
  if (!isInitialised) {
    isInitialised = true;
    const G4String& pname = fParticle->GetParticleName();
    if (pname == "proton") {
      idx = 1;
      xscurrent = xsdata[1];
      fElow = fLowestEnergy;
    } else if (pname == "hydrogen") {
      idx = 0; 
      xscurrent = xsdata[0];
      fElow = fLowestEnergy;
    } else if (pname == "alpha") {
      idx = 1;
      xscurrent = xsdata[2];
      isHelium = true;
      fElow = fLimitEnergy;
    } else if (pname == "alpha+") {
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
    } else if (pname == "helium") {
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
      idx = -1;
      xscurrent = xsdata[1];
      fElow = fLowestEnergy;
    }
    // defined stationary mode
    statCode = G4EmParameters::Instance()->DNAStationary();

    // initialise atomic de-excitation
    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

    if (verbose > 0) {
      G4cout << "### G4DNARuddIonisationExtendedModel::Initialise(..) " << pname 
	     << "/n    idx=" << idx << " isIon=" << isIon
	     << " isHelium=" << isHelium << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::SetParticle(const G4ParticleDefinition* p)
{
  fParticle = p;
  fMass = p->GetPDGMass();
  fMassRate = (isIon) ? CLHEP::proton_mass_c2/fMass : 1.0; 
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

  // ion shoud be stopped - check on kinetic energy and not scaled energy
  if (kinE < fLowestEnergy) { return DBL_MAX; }

  G4double e = kinE*fMassRate;

  G4double sigma = (e > fElow) ? xscurrent->FindValue(e)
    : xscurrent->FindValue(fElow) * e / fElow;

  if (idx == -1) {
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

  G4int shell = SelectShell(kinE*fMassRate);
  G4double bindingEnergy = (useDNAWaterStructure)
    ? waterStructure.IonisationEnergy(shell) : Bj[shell];

  //Si: additional protection if tcs interpolation method is modified
  if (kinE < bindingEnergy) { return; }
  
  G4double esec = SampleElectronEnergy(kinE, shell);
  G4double esum = 0.0;

  // sample deexcitation
  // here we assume that H2O electronic levels are the same as Oxygen.
  // this can be considered true with a rough 10% error in energy on K-shell,
  G4int Z = 8;	
  G4ThreeVector deltaDir = 
    GetAngularDistribution()->SampleDirectionForShell(dpart, esec, Z, shell, couple->GetMaterial());

  // SI: only atomic deexcitation from K shell is considered
  if(fAtomDeexcitation != nullptr && shell == 4) {
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
  auto  dp = new G4DynamicParticle(G4Electron::Electron(), deltaDir, esec);
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
  for (G4int i=0; i<5; ++i) {
    auto ptr = xscurrent->GetComponent(i);
    xs = (e > fElow) ? ptr->FindValue(e) : ptr->FindValue(fElow)*e/fElow;
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

G4double G4DNARuddIonisationExtendedModel::MaxEnergy(G4double kine, G4int shell)
{
  // kinematic limit
  G4double tau = kine/fMass;
  G4double gam = 1.0 + tau;
  G4double emax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0);

  // Initialisation of sampling
  G4double A1, B1, C1, D1, E1, A2, B2, C2, D2;
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
  bEnergy = Bj[shell];
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

G4double G4DNARuddIonisationExtendedModel::SampleElectronEnergy(G4double kine,
                                                                G4int shell)
{
  G4double emax = MaxEnergy(kine, shell);
  // compute cumulative probability function
  G4double step = 1*CLHEP::eV;
  auto nn = (G4int)(emax/step);
  nn = std::min(std::max(nn, 10), 100);
  step = emax/(G4double)nn;

  // find max probability
  G4double pmax = ProbabilityFunction(kine, 0.0, shell);
  //G4cout << "## E(keV)=" << kine/keV << " emax=" << emax/keV
  //       << " pmax(0)=" << pmax << " shell=" << shell << " nn=" << nn << G4endl;

  G4double e0 = 0.0; // energy with max probability
  // 2 areas after point with max probability
  G4double e1 = emax;
  G4double e2 = emax;
  G4double p1 = 0.0;
  G4double p2 = 0.0;
  const G4double f = 0.25;

  // find max probability
  G4double e = 0.0;
  G4double p = 0.0;
  for (G4int i=0; i<nn; ++i) {
    e += step;
    p = ProbabilityFunction(kine, e, shell);
    if (p > pmax) {
      pmax = p;
      e0 = e;
    } else {
      break;
    }
  }
  // increase step to be more effective
  step *= 2.0;
  // 2-nd area
  for (G4int i=0; i<nn; ++i) {
    e += step;
    if (std::abs(e - emax) < step) {
      e1 = emax;
      break;
    }
    p = ProbabilityFunction(kine, e, shell);
    if (p < f*pmax) {
      p1 = p;
      e1 = e;
      break;
    }
  }
  // 3-d area
  if (e < emax) {
    for (G4int i=0; i<nn; ++i) {
      e += step;
      if (std::abs(e - emax) < step) {
        e2 = emax;
	break;
      }
      p = ProbabilityFunction(kine, e, shell);
      if (p < f*p1) {
	p2 = p;
	e2 = e;
        break;
      }
    }
  }
  pmax *= 1.05;
  // regression method with 3 regions
  G4double s0 = pmax*e1;
  G4double s1 = s0 + p1 * (e2 - e1);
  G4double s2 = s1 + p2 * (emax - e2);
  s0 = (s0 == s1) ? 1.0 : s0 / s2;
  s1 = (s1 == s2) ? 1.0 : s1 / s2;

  //G4cout << "pmax=" << pmax << " e1(keV)=" << e1/keV << " p1=" << p1 << " e2(keV)=" << e2/keV
  //	 << " p2=" << p2 << " s0=" << s0 << " s1=" << s1 << " s2=" << s2 << G4endl;

  // sampling
  G4int count = 0;
  G4double ymax, y, deltae;
  for (G4int i = 0; i<100000; ++i) {
    G4double q = G4UniformRand();
    if (q <= s0) {
      ymax = pmax;
      deltae = e1 * q / s0;
    } else if (q <= s1) {
      ymax = p1;
      deltae = e1 + (e2 - e1) * (q - s0) / (s1 - s0);
    } else {
      ymax = p2;
      deltae = e2 + (emax - e2) * (q - s1) / (1.0 - s1);
    }
    y = ProbabilityFunction(kine, deltae, shell);
    //G4cout << "    " << i << ".  deltae=" << deltae/CLHEP::keV 
    //       << " y=" << y << " ymax=" << ymax << G4endl; 
    if (y > ymax && count < 10) {
      ++count;
      G4cout << "G4DNARuddIonisationExtendedModel::SampleElectronEnergy warning: "
	     << fParticle->GetParticleName() << " E(keV)=" << kine/CLHEP::keV
	     << " Edelta(keV)=" << deltae/CLHEP::keV 
	     << " y=" << y << " ymax=" << ymax << " n=" << i << G4endl; 
    }
    if (ymax * G4UniformRand() <= y) {
      return deltae;
    }
  }
  deltae = std::min(e0 + step, 0.5*emax);
  return deltae;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::ProbabilityFunction(G4double kine,
                                                               G4double deltae,
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
  G4double w = deltae/bEnergy;
  G4double x = alphaConst*(w - wc)/v;
  G4double y = (x > -15.) ? 1.0 + G4Exp(x) : 1.0;

  G4double res = CorrectionFactor(kine, shell) * (F1 + w*F2) /
    (fGpow->powN((1. + w)/u, 3) * y);

  if (isHelium) {
    G4double energyTransfer = deltae + bEnergy;
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
  MaxEnergy(e, shell);
  return ProbabilityFunction(e, deltae, shell);
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
