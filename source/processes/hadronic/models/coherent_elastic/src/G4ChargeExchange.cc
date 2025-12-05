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
// G4 Model: Charge and strangness exchange based on G4LightMedia model
//           28 May 2006 V.Ivanchenko
//
// Modified:
// 07-Jun-06 V.Ivanchenko fix problem of rotation of final state
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy, below which S-wave is sampled
// 12-Jun-12 A.Ribon fix warnings of shadowed variables
// 06-Aug-15 A.Ribon migrating to G4Exp, G4Log and G4Pow
//

#include "G4ChargeExchange.hh"
#include "G4ChargeExchangeXS.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4NucleiProperties.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayProducts.hh"
#include "G4NistManager.hh"
#include "G4Fragment.hh"
#include "G4ExcitationHandler.hh"
#include "G4ReactionProductVector.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4HadronicParameters.hh"
#include "G4PhysicsModelCatalog.hh"

namespace
{
  constexpr G4int maxN = 1000;
}

G4ChargeExchange::G4ChargeExchange(G4ChargeExchangeXS* ptr)
  : G4HadronicInteraction("ChargeExchange"),
    fXSection(ptr), fXSWeightFactor(1.0)
{
  lowEnergyLimit = 1.*CLHEP::MeV;
  secID = G4PhysicsModelCatalog::GetModelID( "model_ChargeExchange" );
  nist = G4NistManager::Instance();
  fHandler = new G4ExcitationHandler();
  if (nullptr != fXSection) {
    fXSWeightFactor = 1.0/fXSection->GetCrossSectionFactor();
  }
}

G4ChargeExchange::~G4ChargeExchange()
{
  delete fHandler;
}

G4HadFinalState* G4ChargeExchange::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  auto part = aTrack.GetDefinition();
  G4double ekin = aTrack.GetKineticEnergy();

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();
  
  if (ekin <= lowEnergyLimit) {
    return &theParticleChange;
  }
  theParticleChange.SetWeightChange(fXSWeightFactor);

  G4int projPDG = part->GetPDGEncoding();
  // for hydrogen targets and positive projectile change exchange
  // is not possible on proton, only on deuteron
  if (1 == Z && (211 == projPDG || 321 == projPDG)) { A = 2; } 
  
  if (verboseLevel > 1) {
    G4cout << "G4ChargeExchange for " << part->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z
	   << " A= " << A << " N= " << A - Z
	   << G4endl;
  }

  G4double mass1 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4LorentzVector lv0 = aTrack.Get4Momentum();

  // select final state
  const G4ParticleDefinition* theSecondary =
    fXSection->SampleSecondaryType(part, aTrack.GetMaterial(),
				   Z, A, aTrack.GetTotalEnergy());
  G4int pdg = theSecondary->GetPDGEncoding();

  if (verboseLevel > 1)
    G4cout << "    Secondary " << theSecondary->GetParticleName() << "  pdg=" << pdg << G4endl;

  // omega(782) and f2(1270)
  G4bool isShortLived = (pdg == 223 || pdg == 225);

  // atomic number of the recoil nucleus
  if (projPDG == -211) { --Z; }
  else if (projPDG == 211) { ++Z; }
  else if (projPDG == -321) { --Z; }
  else if (projPDG == 321) { ++Z; }
  else if (projPDG == 130) {
    if (theSecondary->GetPDGCharge() > 0.0) { --Z; }
    else { ++Z; }
  } else {
    // not ready for other projectile
    return &theParticleChange;
  }

  // recoil nucleus
  const G4ParticleDefinition* theRecoil = nullptr;
  if (Z == 0 && A == 1) { theRecoil = G4Neutron::Neutron(); }
  else if (Z == 1 && A == 1) { theRecoil = G4Proton::Proton(); }
  else if (Z == 1 && A == 2) { theRecoil = G4Deuteron::Deuteron(); }
  else if (Z == 1 && A == 3) { theRecoil = G4Triton::Triton(); }
  else if (Z == 2 && A == 3) { theRecoil = G4He3::He3(); }
  else if (Z == 2 && A == 4) { theRecoil = G4Alpha::Alpha(); }

  // check if there is enough energy for the final state
  // and sample mass of produced state
  // sample kinematics
  G4LorentzVector lv1(0.0, 0.0, 0.0, mass1);
  G4LorentzVector lv = lv0 + lv1;
  G4double m0 = lv.mag();
  const G4double mass0 = theSecondary->GetPDGMass();
  G4double mass2 = mass0;
  G4double mass3;
  G4bool ok = false;

  if (verboseLevel > 1) {
    G4cout << " Secondary meson " << theSecondary->GetParticleName()
	   << " mass(MeV)=" << mass2 << " pdg=" << pdg
	   << " Final Z=" << Z << " isShortLived=" << isShortLived
	   << "  " << lv
	   << G4endl;
  }
  // fixed recoil mass
  if (nullptr != theRecoil) {
    mass3 = theRecoil->GetPDGMass();
    ok = (m0 > mass2 + mass3);

    // excited nuclear state
  } else {
    G4double mass30 = G4NucleiProperties::GetNuclearMass(A, Z);
    const G4double eFermi = 10*CLHEP::MeV;
    for (G4int i=0; i<10; ++i) {
      mass3 = mass30 + eFermi*G4UniformRand();
      if (m0 > mass2 + mass3) {
	ok = true;
	break;
      }
    }
  }
  if (isShortLived) {
    const G4double elim = 300*CLHEP::MeV;
    ok = false;
    for (G4int i=0; i<10; ++i) {
      if (SampleMass(mass2, theSecondary->GetPDGWidth(), elim)) {
        if (m0 > mass2 + mass3) {
	  ok = true;
	  break;
	}
      }
    }
  }

  // not possible kinematically
  if (!ok) { return &theParticleChange; }

  G4double e2 = (m0*m0 + mass2*mass2 - mass3*mass3)/(2*m0);
  G4double momentumCMS = std::sqrt(e2*e2 - mass2*mass2);
  G4double tmax = 4*momentumCMS*momentumCMS; 

  // for projectile pion t depends on final state
  G4double t;
  if (fXSection->isPion()) {
    t = fXSection->SampleTforPion(aTrack.GetTotalEnergy(), tmax);
  }
  else {
    t = SampleT(theSecondary, A, tmax);
  } 
   
  G4double phi = G4UniformRand()*CLHEP::twopi;
  G4double cost = 1. - 2.0*t/tmax;

  // if cos(theta) negative, there is a numerical problem
  // instead of making scattering backward, make in this case
  // no scattering
  if (std::abs(cost) > 1.0) { cost = 1.0; }

  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));

  if (verboseLevel > 1) {
    G4cout << " t= " << t << " tmax(GeV^2)= " << tmax/(GeV*GeV) 
	   << " cos(t)=" << cost << " sin(t)=" << sint << G4endl;
  }
  G4LorentzVector lv2(momentumCMS*sint*std::cos(phi),
		      momentumCMS*sint*std::sin(phi),
		      momentumCMS*cost, e2);

  // kinematics in the final state, may be a warning should be added if 
  G4ThreeVector bst = lv.boostVector();
  lv2.boost(bst);
  lv -= lv2;
  if (lv.e() < mass3) {
    lv.setE(mass3);
  }

  // prepare secondary particles
  theParticleChange.SetStatusChange(stopAndKill);
  theParticleChange.SetEnergyChange(0.0);
  theParticleChange.SetWeightChange(fXSWeightFactor);

  if (!isShortLived) {
    auto aSec = new G4DynamicParticle(theSecondary, lv2);
    theParticleChange.AddSecondary(aSec, secID);
  } else {
    auto channel = theSecondary->GetDecayTable()->SelectADecayChannel(mass2);
    auto products = channel->DecayIt(mass2);
    G4ThreeVector bst1 = lv2.boostVector();
    G4int N = products->entries();
    for (G4int i=0; i<N; ++i) {
      auto p = (*products)[i];
      auto lvp = p->Get4Momentum();
      lvp.boost(bst1);
      auto pnew = new G4DynamicParticle(*p);
      pnew->Set4Momentum(lvp);
      theParticleChange.AddSecondary(pnew, secID);
    }
    delete products;
  }

  // recoil is a stable isotope
  if (nullptr != theRecoil) {
    auto aRec = new G4DynamicParticle(theRecoil, lv);
    theParticleChange.AddSecondary(aRec, secID);
  } else {
    // recoil is a fragment, which may be unstable
    G4Fragment frag(A, Z, lv);
    auto products = fHandler->BreakItUp(frag);
    for (auto & prod : *products) {
      auto dp = new G4DynamicParticle(prod->GetDefinition(), prod->GetMomentum());
      theParticleChange.AddSecondary(dp, secID);
      delete prod;
    }
    delete products;
  }
  return &theParticleChange;
}

G4double G4ChargeExchange::SampleT(const G4ParticleDefinition*,
                                   const G4int A, const G4double ltmax) const
{
  const G4double GeV2 = CLHEP::GeV*CLHEP::GeV;
  const G4double numLimit = 18.;

  G4double tmax = ltmax/GeV2;
  if (verboseLevel > 1) {
    G4cout << "G4ChargeExchange::SampleT tmax(GeV^2)=" << tmax << G4endl;
  }

  G4double aa, bb, cc, dd;
  G4Pow* g4pow = G4Pow::GetInstance();
  G4double a13 = g4pow->Z13(A);
  if (A <= 62) {
    aa = (A*A);
    bb = 14.5*a13*a13;
    cc = 1.4*a13;
    dd = 10.;
  } else {
    aa = g4pow->powZ(A, 1.33);
    bb = 60.*a13;
    cc = 0.4*g4pow->powZ(A, 0.40);
    dd = 10.;
  }
  G4double q1 = 1.0 - G4Exp(-std::min(bb*tmax, numLimit));
  G4double q2 = 1.0 - G4Exp(-std::min(dd*tmax, numLimit));
  G4double s1 = q1*aa;
  G4double s2 = q2*cc;
  if ((s1 + s2)*G4UniformRand() < s2) {
    q1 = q2;
    bb = dd;
  }
  return -GeV2*G4Log(1.0 - G4UniformRand()*q1)/bb;
}

G4bool G4ChargeExchange::SampleMass(G4double& M, const G4double G,
				    const G4double elim)
{
  // +- 4 width but above 2 pion mass
  G4double e1 = std::max(M - 4*G, elim);
  G4double e2 = M + 4*G - e1;
  if (e2 <= 0.0) { return false; }
  G4double M2 = M*M;
  G4double MG2 = M2*G*G;

  // sampling Breit-Wigner function
  for (G4int i=0; i<maxN; ++i) {
    G4double e = e1 + e2*G4UniformRand();
    G4double x = e*e - M2;
    G4double y = MG2/(x*x + MG2);
    if (y >= G4UniformRand()) {
      M = e;
      return true;
    }
  }
  return false;
}
