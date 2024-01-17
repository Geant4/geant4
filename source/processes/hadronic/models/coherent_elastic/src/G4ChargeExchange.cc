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

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4HadronicParameters.hh"
#include "G4PhysicsModelCatalog.hh"


G4ChargeExchange::G4ChargeExchange(G4ChargeExchangeXS* ptr)
  : G4HadronicInteraction("ChargeExchange"),
    fXSection(ptr)
{
  lowEnergyLimit = 1.*CLHEP::MeV;
  secID = G4PhysicsModelCatalog::GetModelID( "model_ChargeExchange" );
}

G4HadFinalState* G4ChargeExchange::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  auto part = aTrack.GetDefinition();
  G4int pdg = part->GetPDGEncoding();
  G4double ekin = aTrack.GetKineticEnergy();

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  if (ekin <= lowEnergyLimit) {
    return &theParticleChange;
  }

  G4int projPDG = part->GetPDGEncoding();
  if (verboseLevel > 1)
    G4cout << "G4ChargeExchange for " << part->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z
	   << " A= " << A << " N= " << A - Z
	   << G4endl;

  G4double mass1 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4LorentzVector lv0 = aTrack.Get4Momentum();
  G4double etot = mass1 + lv0.e();

  // select final state
  const G4ParticleDefinition* theRecoil = nullptr;
  const G4ParticleDefinition* theSecondary =
    fXSection->SampleSecondaryType(part, Z, A);

  // atomic number of the recoil nucleus
  if (pdg == -211) { --Z; }
  else if (pdg == 211) { ++Z; }
  else if (pdg == -321) { --Z; }
  else if (pdg == 321) { ++Z; }
  else if (pdg == 130) {
    if (theSecondary->GetPDGCharge() > 0.0) { --Z; }
    else { ++Z; }
  } else {
    // not ready for other projectile
    return &theParticleChange;
  }

  if (Z == 0 && A == 1) theRecoil = G4Neutron::Neutron();
  else if (Z == 1 && A == 1) theRecoil = G4Proton::Proton();
  else if (Z == 1 && A == 2) theRecoil = G4Deuteron::Deuteron();
  else if (Z == 1 && A == 3) theRecoil = G4Triton::Triton();
  else if (Z == 2 && A == 3) theRecoil = G4He3::He3();
  else if (Z == 2 && A == 4) theRecoil = G4Alpha::Alpha();
  else {
    theRecoil = G4ParticleTable::GetParticleTable()
      ->GetIonTable()->GetIon(Z, A, 0.0);
  }
  if (nullptr == theRecoil) { return &theParticleChange; }

  G4double mass2 = theSecondary->GetPDGMass();
  G4double mass3 = theRecoil->GetPDGMass();

  if (etot <= mass2 + mass3) {
    // not possible kinematically
    return &theParticleChange;
  }

  // sample kinematics
  G4LorentzVector lv1(0.0, 0.0, 0.0, mass1);
  G4LorentzVector lv = lv0 + lv1;
  G4ThreeVector bst = lv.boostVector();
  G4double ss = lv.mag2();

  // tmax = 4*momCMS^2
  G4double e2 = ss + mass2*mass2 - mass3*mass3;
  G4double tmax = e2*e2/ss - 4*mass2*mass2;

  G4double t = SampleT(theSecondary, A, tmax);
  
  G4double phi  = G4UniformRand()*CLHEP::twopi;
  G4double cost = 1. - 2.0*t/tmax;

  if (cost > 1.0) { cost = 1.0; }
  else if(cost < -1.0) { cost = -1.0; } 

  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  if (verboseLevel>1) {
    G4cout << " t= " << t << " tmax(GeV^2)= " << tmax/(GeV*GeV) 
	   << " cos(t)=" << cost << " sin(t)=" << sint << G4endl;
  }
  G4double momentumCMS = 0.5*std::sqrt(tmax);
  G4LorentzVector lv2(momentumCMS*sint*std::cos(phi),
		      momentumCMS*sint*std::sin(phi),
		      momentumCMS*cost,
		      std::sqrt(momentumCMS*momentumCMS + mass2*mass2));

  // kinematics in the final stae, may be a warning should be added if 
  lv2.boost(bst);
  if (lv2.e() < mass2) {
    lv2.setE(mass2);
  }
  lv -= lv2;
  if (lv.e() < mass3) {
    lv.setE(mass3);
  }

  // prepare secondary particles
  theParticleChange.SetStatusChange(stopAndKill);
  theParticleChange.SetEnergyChange(0.0);

  G4DynamicParticle * aSec = new G4DynamicParticle(theSecondary, lv2);
  theParticleChange.AddSecondary(aSec, secID);

  G4DynamicParticle * aRec = new G4DynamicParticle(theRecoil, lv);
  theParticleChange.AddSecondary(aRec, secID);
  return &theParticleChange;
}

G4double G4ChargeExchange::SampleT(const G4ParticleDefinition*,
                                   G4int A, G4double tmax) const
{
  G4double aa, bb, cc, dd;
  G4Pow* g4pow = G4Pow::GetInstance();
  if (A <= 62.) {
    aa = g4pow->powZ(A, 1.63);
    bb = 14.5*g4pow->powZ(A, 0.66);
    cc = 1.4*g4pow->powZ(A, 0.33);
    dd = 10.;
  } else {
    aa = g4pow->powZ(A, 1.33);
    bb = 60.*g4pow->powZ(A, 0.33);
    cc = 0.4*g4pow->powZ(A, 0.40);
    dd = 10.;
  }
  G4double x1 = (1.0 - G4Exp(-tmax*bb))*aa/bb;
  G4double x2 = (1.0 - G4Exp(-tmax*dd))*cc/dd;
  
  G4double t;
  G4double y = bb;
  if(G4UniformRand()*(x1 + x2) < x2) y = dd;

  const G4int maxN = 10000;
  G4int count = 0;
  do {
    t = -G4Log(G4UniformRand())/y;
  } while ( (t > tmax) && ++count < maxN );
  /* Loop checking, 10.08.2015, A.Ribon */
  if ( count >= maxN ) {
    t = 0.0;
  }
  return t;
}
