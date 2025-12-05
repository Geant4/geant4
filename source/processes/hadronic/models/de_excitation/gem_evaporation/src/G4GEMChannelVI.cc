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
// GEM de-excitation model
// by V. Ivanchenko (July 2019)
//

#include "G4GEMChannelVI.hh"
#include "G4GEMProbabilityVI.hh"
#include "G4VCoulombBarrier.hh"
#include "G4CoulombBarrier.hh"
#include "G4DeexPrecoUtility.hh"
#include "G4PairingCorrection.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4NucleiProperties.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicsModelCatalog.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4InterfaceToXS.hh"
#include "G4IsotopeList.hh"
#include "G4NuclearRadii.hh"

namespace
{
  const G4double dExc = 1.0*CLHEP::MeV; 
  const G4double limE = 1.0*CLHEP::MeV; // low-energy limit for neutrons
  const G4double kmin = 20*CLHEP::keV;  // low-energy limit on primary kinetic energy
  const G4int nProbMax = 10;
  G4double prob[nProbMax] = {0.0};  
}

G4GEMChannelVI::G4GEMChannelVI(G4int theA, G4int theZ)
  : evapA(theA), evapZ(theZ)
{
  nData = G4NuclearLevelData::GetInstance();
  pairingCorrection = nData->GetPairingCorrection();
  if (evapZ > 2) { lManagerEvap = nData->GetLevelManager(evapZ, evapA); }
  fEvapMass = G4NucleiProperties::GetNuclearMass(evapA, evapZ);
  fEvapMass2 = fEvapMass*fEvapMass;

  cBarrier = new G4CoulombBarrier(evapA, evapZ);

  fTolerance = 10*CLHEP::keV;
  fCoeff = fEvapMass/((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc));

  std::ostringstream ss;
  ss << "GEMVI_" << "Z" << evapZ << "_A" << evapA;
  fModelName = ss.str();
  
  fNeutron = G4Neutron::Neutron();
  fProton = G4Proton::Proton();

  secID = G4PhysicsModelCatalog::GetModelID("model_G4GEMChannelVI");
  const G4ParticleDefinition* part = nullptr;
  if (evapZ == 0 && evapA == 1) {
    indexC = 0;
    fCoeff *= 2.0;
    part = fNeutron;
  } else if (evapZ == 1 && evapA == 1) {
    indexC = 1;
    fCoeff *= 2.0;
    part = fProton;
  } else if (evapZ == 1 && evapA == 2) {
    indexC = 2;
    fCoeff *= 3.0;
    part = G4Deuteron::Deuteron();
  } else if (evapZ == 1 && evapA == 3) {
    indexC = 3;
    fCoeff *= 2.0;
    part = G4Triton::Triton();
  } else if (evapZ == 2 && evapA == 3) {
    indexC = 4;
    fCoeff *= 2.0;
    part = G4He3::He3();
  } else if (evapZ == 2 && evapA == 4) {
    indexC = 5;
    part = G4Alpha::Alpha();
  } else {
    G4int N = evapA - evapZ;
    fCoeff *= (1 + (evapZ - 2*(evapZ/2)))*(1 + (N - 2*(N/2)));
  }
  g4pow = G4Pow::GetInstance();
  
  G4double de = (0 == indexC) ? 0.15*CLHEP::MeV : 0.25*CLHEP::MeV;
  InitialiseIntegrator(0.01, 0.25, 1.3, de, 0.1*CLHEP::MeV, 2*CLHEP::MeV);

  if (indexC <= 6) { fXSection = new G4InterfaceToXS(part, indexC); }
}

G4GEMChannelVI::~G4GEMChannelVI()
{
  delete cBarrier;
  delete fXSection;
}

void G4GEMChannelVI::Initialise()
{
  G4VEvaporationChannel::Initialise(); 
}

G4double G4GEMChannelVI::GetEmissionProbability(G4Fragment* fragment)
{
  fragZ = fragment->GetZ_asInt();
  fragA = fragment->GetA_asInt();
  resZ = fragZ - evapZ;
  resA = fragA - evapA;
  // G4cout << "G4GEMChannelVI::GetEmissionProbability Z=" << evapZ <<  " A=" << evapA << " resZ=" << resZ << " resA=" << resA << G4endl;
  // to avoid double counting
  if (resA < evapA || resA < resZ || resZ < 1 ||
      (resA == evapA && resZ < evapZ)) { return 0.0; }

  fFragExc = fragment->GetExcitationEnergy();
  fMass = fragment->GetGroundStateMass() + fFragExc;
  fResMass = G4NucleiProperties::GetNuclearMass(resA, resZ);
  fResA13 = g4pow->Z13(resA);
  xsfactor = g4pow->Z23(fragA)/g4pow->Z23(resA);

  // limit for the case when both evaporation and residual 
  // fragments are in ground states
  if (fMass <= fEvapMass + fResMass) { return 0.0; } 

  a0 = G4DeexPrecoUtility::LevelDensity(fragZ, fragA, indexC);
  delta0 = nData->GetPairingCorrection(fragZ, fragA);
  delta1 = nData->GetPairingCorrection(resZ, resA);
  fE0 = std::max(fFragExc - delta0, 0.0);
  
  if (indexC > 0) {
    bCoulomb = cBarrier->GetCoulombBarrier(resA, resZ, fFragExc);
  }
  G4double elim = 0.5*bCoulomb;
  G4double de = fMass - fEvapMass - fResMass - elim;
  if (de < fTolerance) { return 0.0; }
  nProbEvap = 1;
  fDeltaEvap = de;
  if (7 == indexC) {
    G4int n = (G4int)(de/dExc) + 1;
    nProbEvap = std::min(n, nProbMax);
    if (nProbEvap > 1) { fDeltaEvap /= (G4double)(nProbEvap - 1); }
  }

  if (2 < fVerbose) {
    G4cout << "## G4GEMChannelVI::GetEmissionProbability fragZ="
	   << fragZ << " fragA=" << fragA << " Z=" << evapZ << " A=" << evapA
	   << " Eex(MeV)=" << fFragExc << " nProbEvap=" << nProbEvap
	   << " nProbRes=" << nProbRes << " CB=" << bCoulomb
	   << " Elim=" << fEnergyLimitXS << " XSfac=" << xsfactor << G4endl;
  }

  // m1 is the mass of emitted excited fragment
  // e2 - free energy in the 2-body decay
  G4double sump = 0.0;
  for (G4int i = 0; i < nProbEvap; ++i) {
    fEvapExc = fDeltaEvap*i;
    G4double m1 = fEvapMass + fEvapExc;
    G4double e2 = fMass - m1 - fResMass;
    e2 = std::max(e2, 0.0);
    if (e2 <= elim + fTolerance) {
      nProbEvap = i + 1;
      prob[i] = sump;
      break;
    }
    G4double p = ComputeIntegral(elim, e2);
    sump += p;
    prob[i] = sump;
    if (2 < fVerbose) {
      G4cout << i << ". e1=" << elim << " e2=" << e2 << " e2-e1="
	     << e2 - elim << " fEvapExc=" << fEvapExc
	     << " Probability=" << p << G4endl;
    }
  }
  if (nProbEvap > 1) { sump /= (G4double)nProbEvap; }
  return sump;
}

G4double G4GEMChannelVI::ProbabilityDensityFunction(G4double e)
{
  // e is free energy 
  G4double m1 = fEvapMass + fEvapExc;
  fResExc = fMass - m1 - fResMass - e;
  if (fResExc < 0.0 || 0.0 == e) { return 0.0; }
  fE1 = std::max(fResExc - delta1, 0.0);
  a1 = G4DeexPrecoUtility::LevelDensity(resZ, resA, indexC);
  G4double m2 = fResMass + fResExc;
  G4double elab = 0.5*(fMass + m1 + m2)*(fMass - m1 - m2)/m2;
  G4double xs = CrossSection(elab);
  G4double res =
    fCoeff*G4Exp(2.0*(std::sqrt(a1*fE1) - std::sqrt(a0*fE0)))*e*xs;

  //G4cout << "e=" << e << " elab=" << elab << " xs(mb)="
  //	 << xs/CLHEP::millibarn << " prob=" << res << G4endl;
  return res;
}

G4double G4GEMChannelVI::CrossSection(G4double e)
{
  G4int Z = std::min(resZ, ZMAXNUCLEARDATA);
  G4double corr;
  G4double e1 = std::max(e, kmin);
  if (e1 < 0.5*bCoulomb) {
    recentXS = 0.0;
    return recentXS;
  }
  if (indexC <= 5) {
    G4double e2 = 2*bCoulomb;
    if (0 == indexC) {
      e2 = lowEnergyLimitMeV[Z];
      if (e2 == 0.0) { e2 = limE; }
    }
    e1 = std::max(e1, e2);
    corr = G4DeexPrecoUtility::CorrectionFactor(indexC, evapZ, fResA13, bCoulomb, e);
    recentXS = fXSection->GetElementCrossSection(e1, Z)/CLHEP::millibarn;
  } else {
    G4double tR = G4NuclearRadii::Radius(resZ, resA);  
    G4double pR = G4NuclearRadii::Radius(evapZ, evapA);
    corr = G4DeexPrecoUtility::CorrectionFactor(indexC, Z, fResA13, bCoulomb, e);

    // geometrical x-section
    recentXS = CLHEP::pi*(pR + tR)*(pR + tR);
  }
  recentXS *= corr;
  return recentXS;
}

G4Fragment* G4GEMChannelVI::EmittedFragment(G4Fragment* theNucleus)
{
  // assumed, that TotalProbability(...) was already called
  // if value iz zero no possiblity to sample final state
  G4Fragment* evFragment = nullptr;
  lManagerRes = nData->GetLevelManager(resZ, resA);
  G4double e2 = fMass - fEvapMass - fResMass;
  fEvapExc = 0.0;

  // sample excitation of the evaporation fragment
  if (nProbEvap > 1) {
    G4double q = prob[nProbEvap - 1];
    if (q > 0.0) {
      q *= G4UniformRand();
      for (G4int i=0; i < nProbEvap; ++i) {
	if (q <= prob[i]) {
	  if (0 == i) { break; }
	  G4double e1 = fDeltaEvap*((i - 1) + (q - prob[i - 1])/(prob[i] - prob[i - 1]));
	  fEvapExc = CorrectExcitation(e1, lManagerEvap);
	  e2 -= fEvapExc;
	  e2 = std::max(e2, 0.0);
	}
      }
    }
  }
  if (ComputeIntegral(0.5*bCoulomb, e2) <= 0.0) { return evFragment; }

  // sample free energy
  G4double e = SampleValue();
  // compute excitation of the residual fragment
  fResExc = CorrectExcitation(e2 - e, lManagerRes);

  // final kinematics
  G4double m1 = fEvapMass + fEvapExc;
  G4double m2 = fResMass + fResExc;

  G4double ekin = 0.5*e*(e + 2*m2)/(e + m1 + m2);
  G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*m1))
		     *G4RandomDirection(), ekin + m1);
  G4LorentzVector lv0 = theNucleus->GetMomentum();
  lv.boost(lv0.boostVector());
  evFragment = new G4Fragment(evapA, evapZ, lv);
  evFragment->SetCreatorModelID(secID);

  // residual
  lv0 -= lv;
  theNucleus->SetZandA_asInt(resZ, resA);
  theNucleus->SetMomentum(lv0);
  theNucleus->SetCreatorModelID(secID);
  
  return evFragment;  
}

G4double
G4GEMChannelVI::CorrectExcitation(G4double exc, const G4LevelManager* man)
{
  if (exc <= 0.0 || nullptr == man) { return 0.0; }
  std::size_t idx = man->NearestLevelIndex(exc);

  // choose ground state
  if (0 == idx) { return 0.0; }

  // possible discrete level 
  G4double elevel = man->LevelEnergy(idx);
  std::size_t ntrans{0};
  if (std::abs(elevel - exc) < fTolerance) {
    auto level = man->GetLevel(idx);
    if (nullptr != level) { 
      ntrans = level->NumberOfTransitions();
      G4int idxfl = man->FloatingLevel(idx); 
      // for floating level check levels with the same energy
      if (idxfl > 0) {
	auto newlevel = man->GetLevel(idx - 1);
	G4double newenergy = man->LevelEnergy(idx - 1);
	if (nullptr != newlevel && std::abs(elevel - newenergy) < fTolerance) {
	  std::size_t newntrans = newlevel->NumberOfTransitions();
	  if (newntrans > 0) {
	    elevel = newenergy;
	    ntrans = newntrans;
          }
        }
      }
      if (0 < ntrans) { return elevel; }
    }
  }
  return exc;
}

const G4String& G4GEMChannelVI::ModelName() const
{
  return fModelName;
}

void G4GEMChannelVI::Dump() const
{}



