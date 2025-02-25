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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4ChargeExchangeXS
//

#include "G4ChargeExchangeXS.hh"
#include "G4DynamicParticle.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4IsotopeList.hh"
#include "G4HadronicParameters.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"  
#include "G4Pow.hh"

#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4ParticleTable.hh"

namespace {
  // V. Lyubovitsky parameterisation
  const G4double piA[5] = {430., 36.,  1.37,  2.0,  60.};  // A
  const G4double pAP[5] = {1.04, 1.26, 1.35,  0.94, 0.94}; // 2 - 2alphaP
  const G4double pC0[5] = {12.7, 6.0,  6.84,  6.5,  8.0};  // c0
  const G4double pC1[5] = {1.57, 1.6,  1.7,   1.23, 2.6};  // c1
  const G4double pG0[5] = {2.55, 4.6,  3.7,   5.5,  4.6};  // g0
  const G4double pG1[5] = {-0.23, -0.5,  0.,    0., -2.};  // g1

  const G4double beta_prime_pi = 0.0410;
}

G4ChargeExchangeXS::G4ChargeExchangeXS() 
{
  if (verboseLevel > 1) {
    G4cout  << "G4ChargeExchangeXS::G4ChargeExchangeXS" << G4endl;
  }
  g4calc = G4Pow::GetInstance();
  auto table = G4ParticleTable::GetParticleTable();
  const G4String nam[5] = {"pi0", "eta", "eta_prime", "omega", "f2(1270)"};
  for (G4int i=0; i<5; ++i) {
    fPionSecPD[i] = table->FindParticle(nam[i]);
    if (nullptr == fPionSecPD[i]) {
      G4ExceptionDescription ed;
      ed << "### meson " << nam[i] << " is not found out in the particle table";
      G4Exception("G4ChargeExchangeXS::G4ChargeExchangeXS()","had044",
                  FatalException, ed,"");
    }
  }
}

// Print the information of this .cc file
void G4ChargeExchangeXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ChargeExchangeXS calculates charge exchange cross section for "
          << "pi+, pi-, K+, K-, KL\n";
}

G4bool G4ChargeExchangeXS::IsElementApplicable(const G4DynamicParticle*,
                                               G4int, const G4Material*)
{
  return true;
}

G4double
G4ChargeExchangeXS::GetElementCrossSection(const G4DynamicParticle* aParticle, 
				           G4int ZZ, const G4Material* mat)  
{
  G4double result = 0.0;
  const G4double pE = aParticle->GetTotalEnergy();
  if (pE <= fEnergyLimit) { return result; }
  
  auto part = aParticle->GetDefinition();
  G4int pdg = part->GetPDGEncoding();   

  // Get or calculate the proton mass, particle mass, and s(Lorentz invariant) 
  G4double tM = CLHEP::proton_mass_c2;
  G4double pM = part->GetPDGMass();
  G4double lorentz_s = tM*tM + 2*tM*pE +  pM*pM;
  if (lorentz_s <= (tM + pM)*(tM + pM)) { return result; }

  const G4int Z = std::min(ZZ, ZMAXNUCLEARDATA);
  const G4int A = G4lrint(aeff[Z]);

  if (verboseLevel > 1) {
    G4cout << "### G4ChargeExchangeXS: " << part->GetParticleName()
	   << " Z=" << Z << " A=" << A << " Etot(GeV)=" << pE/CLHEP::GeV
	   << " s(GeV^2)=" << lorentz_s/(CLHEP::GeV*CLHEP::GeV) << G4endl;
  }

  // For unit conversion 
  const G4double inv1e7 = 0.1/(CLHEP::GeV*CLHEP::GeV);
  const G4double fact = 1e-30*CLHEP::cm2;
  const G4double pfact = 0.1/CLHEP::GeV;
  const G4double kfact = 56.3*fact;
  const G4double csmax = 1e-16;

  // The approximation of Glauber-Gribov formula -> extend it from interaction with 
  // proton to nuclei Z^(2/3). The factor g4calc->powA(A,-beta_prime_pi*G4Log(A))
  // takes into account absorption of pi0 and eta 

  // pi- + p -> n + meson (0- pi0, 1- eta, 2- eta', 3- omega, 4- f2(1270))
  if (pdg == -211) {
    G4double z23 = g4calc->Z23(Z);
    G4double x = lorentz_s*inv1e7;
    G4double logX = G4Log(x);
    G4double logA = g4calc->logZ(A);
    G4double xf = g4calc->powZ(A, -beta_prime_pi*logA);
    G4double sum = 0.0;
    for (G4int i=0; i<5; ++i) {
      G4double xg = std::max(1.0 + pG0[i] + pG1[i]*logX, 0.0);
      G4double xc = std::max(pC0[i] + pC1[i]*logX, csmax);
      G4double xs = z23*piA[i]*g4calc->powA(x, -pAP[i])*xf*xg/xc;
      sum += xs;
      fXSecPion[i] = sum;
    }
    result = sum*fact;
  }

  // pi+ + n -> p + meson (0- pi0, 1- eta, 2- eta', 3- omega, 4- f2(1270))
  else if (pdg == 211) {
    G4double n23 = g4calc->Z23(A - Z);
    G4double x = lorentz_s*inv1e7;
    G4double logX = G4Log(x);
    G4double logA = g4calc->logZ(A);
    G4double xf = g4calc->powZ(A, -beta_prime_pi*logA);

    // hydrogen target case Z = A = 1
    // the cross section is defined by fraction of deuteron and tritium
    if (1 == Z) { n23 = ComputeDeuteronFraction(mat); }
    G4double sum = 0.0;
    for (G4int i=0; i<5; ++i) {
      G4double xg = std::max(1.0 + pG0[i] + pG1[i]*logX, 0.0);
      G4double xc = std::max(pC0[i] + pC1[i]*logX, csmax);
      G4double xs = n23*piA[i]*g4calc->powA(x, -pAP[i])*xf*xg/xc;
      sum += xs;
      fXSecPion[i] = sum;
    }
    result = sum*fact;
  }

  // Kaon x-sections depend on the primary particles momentum
  // K- + p -> Kbar + n
  else if (pdg == -321) {
    G4double p_momentum = std::sqrt(pE*pE - pM*pM)*pfact;
    result = g4calc->Z23(Z)*g4calc->powA(p_momentum, -1.60)*kfact;
  }

  // K+ + n -> Kbar + p
  else if (pdg == 321) {
    G4double p_momentum = std::sqrt(pE*pE - pM*pM)*pfact;
    G4double n23 = g4calc->Z23(A-Z);   
    // hydrogen target case Z = A = 1
    // the cross section is defined by fraction of deuteron and tritium
    if (1 == Z) { n23 = ComputeDeuteronFraction(mat); }
    result = n23*g4calc->powA(p_momentum, -1.60)*kfact;
  }

  // KL 
  else if (pdg == 130) {
    // Cross section of KL = 0.5*(Cross section of K+ + Cross section of K-)
    const G4double p_momentum = std::sqrt(pE*pE - pM*pM)*pfact;
    result = 0.5*(g4calc->Z23(Z) + g4calc->Z23(A-Z))*
      g4calc->powA(p_momentum, -1.60)*kfact;
  }
  result *= fFactor;
  if (verboseLevel > 1) {
    G4cout  << "   Done for " << part->GetParticleName() << " Etot(GeV)=" << pE/CLHEP::GeV
	    << " res(mb)=" << result/CLHEP::millibarn << G4endl;
  }
  return result;
}

const G4ParticleDefinition*
G4ChargeExchangeXS::SampleSecondaryType(const G4ParticleDefinition* part,
                                        const G4int Z, const G4int A)
{
  const G4ParticleDefinition* pd = nullptr;
  G4int pdg = std::abs(part->GetPDGEncoding());  

  // pi- + p /  pi+ + n  
  if (pdg == 211) {
    pd = fPionSecPD[0];
    G4double x = fXSecPion[4]*G4UniformRand();
    for (G4int i=0; i<5; ++i) {
      if (x <= fXSecPion[i]) {
        pd = fPionSecPD[i];
	break;
      }
    }
  }

  // K- + p /  K+ + n 
  // Equal opportunity of producing k-short and k-long
  else if (pdg == 321) {
    if (G4UniformRand() >= 0.5) {
      pd = G4KaonZeroLong::KaonZeroLong();
    }
    else {
      pd = G4KaonZeroShort::KaonZeroShort();
    }
  }

  // KL + nucleus
  else if (pdg == 130) {
    G4double prob = (G4double)Z/(G4double)A;
    if (G4UniformRand() >= prob) {
      pd = G4KaonMinus::KaonMinus();
    }
    else {
      pd = G4KaonPlus::KaonPlus();
    }
  }

  return pd;
} 

G4double
G4ChargeExchangeXS::ComputeDeuteronFraction(const G4Material* mat)
{
  for (auto const & elm : *mat->GetElementVector()) {
    if (1 == elm->GetZasInt()) {
      G4double ab = 0.0;
      const G4int nIso = (G4int)elm->GetNumberOfIsotopes();
      const G4double* abu = elm->GetRelativeAbundanceVector();
      for (G4int j = 0; j < nIso; ++j) {
	auto const iso = elm->GetIsotope(j);
	ab += (iso->GetN() - iso->GetZ())*abu[j];
      }
      return ab;
    }
  }
  return 0.0;
}
