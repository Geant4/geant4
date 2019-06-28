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
// Geant4 Header : G4HadronElastic
//
// Author : V.Ivanchenko 29 June 2009 (redesign old elastic model)
//  

#include "G4HadronElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4HadronicParameters.hh"


G4HadronElastic::G4HadronElastic(const G4String& name) 
  : G4HadronicInteraction(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  lowestEnergyLimit= 1.e-6*eV;  

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();
}

G4HadronElastic::~G4HadronElastic()
{}


void G4HadronElastic::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4HadronElastic is the base class for all hadron-nucleus\n" 
          << "elastic scattering models except HP.\n" 
          << "By default it uses the Gheisha two-exponential momentum\n"
	  << "transfer parameterization.  The model is fully relativistic\n"
	  << "as opposed to the original Gheisha model which was not.\n"
	  << "This model may be used for all long-lived hadrons at all\n"
	  << "incident energies but fit the data only for relativistic scattering.\n";
}

G4HadFinalState* G4HadronElastic::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();

  // no scattering below the limit
  if(ekin <= lowestEnergyLimit) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(0.,0.,1.);
    return &theParticleChange;
  }

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  G4double plab = std::sqrt(ekin*(ekin + 2.0*m1));

  if (verboseLevel>1) {
    G4cout << "G4HadronElastic: " 
	   << aParticle->GetDefinition()->GetParticleName() 
	   << " Plab(GeV/c)= " << plab/GeV  
	   << " Ekin(MeV) = " << ekin/MeV 
	   << " scattered off Z= " << Z 
	   << " A= " << A 
	   << G4endl;
  }

  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double e1 = m1 + ekin;
  G4LorentzVector lv(0.0,0.0,plab,e1+mass2);
  G4ThreeVector bst = lv.boostVector();
  G4double momentumCMS = plab*mass2/std::sqrt(m1*m1 + mass2*mass2 + 2.*mass2*e1);

  pLocalTmax = 4.0*momentumCMS*momentumCMS;

  // Sampling in CM system
  G4double t    = SampleInvariantT(theParticle, plab, Z, A);
  G4double phi  = G4UniformRand()*CLHEP::twopi;
  G4double cost = 1. - 2.0*t/pLocalTmax;

  // For the very rare cases where cos(theta) is greater than 1 or smaller than -1,
  // print some debugging information via a "JustWarning" exception, and safely
  // return (simply setting "cost=1.0" or "cost=-1.0" can sometimes cause a crash,
  // due to numerical imprecisions, e.g. 3-momentum = (0.0, 0.0, 0.0) but
  // Ekin very small but not 0.0).
  if ( std::abs( cost ) > 1.0 ) {
    G4ExceptionDescription ed;
    ed << " LARGE cost ! cost=" << cost << " for " << aParticle->GetDefinition()->GetParticleName() 
       << " ekin=" << ekin << " MeV" << "  on (Z,A)=(" << Z << "," << A << ")" << G4endl;
    if ( cost > 1.0 ) { 
      // We assume here no interaction and let the projectile keep going unchanged.
      theParticleChange.SetEnergyChange( ekin );
      theParticleChange.SetMomentumChange( aParticle->Get4Momentum().vect().unit() );
      ed << "\t No interaction: the projectile keeps going unchanged!" << G4endl;
      G4Exception( "G4HadronElastic::ApplyYourself", "hadEla001", JustWarning, ed );
      return &theParticleChange;
    } else {  // cost < -1.0 ) { 
      // We assume here that the projectile stops and its energy is deposited locally
      // (for simplicity, given that this condition should happen rarely, we neglect
      // the recoil of the target nucleus).
      theParticleChange.SetEnergyChange( 0.0 );
      theParticleChange.SetLocalEnergyDeposit( ekin );
      ed << "\t Projectile stops and its energy is deposited locally:" << G4endl
         << "\t neglected recoil of the target nucleus!" << G4endl;
      G4Exception( "G4HadronElastic::ApplyYourself", "hadEla002", JustWarning, ed );
      return &theParticleChange;
    }
  }

  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  if (verboseLevel>1) {
    G4cout << " t= " << t << " tmax(GeV^2)= " << pLocalTmax/(GeV*GeV) 
	   << " Pcms(GeV)= " << momentumCMS/GeV << " cos(t)=" << cost 
	   << " sin(t)=" << sint << G4endl;
  }
  G4LorentzVector nlv1(momentumCMS*sint*std::cos(phi),
		       momentumCMS*sint*std::sin(phi),
                       momentumCMS*cost,
		       std::sqrt(momentumCMS*momentumCMS + m1*m1));

  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;
  if (verboseLevel > 1) {
    G4cout <<"G4HadronElastic: m= " << m1 << " Efin(MeV)= " << eFinal 
	   << " 4-M Final: " << nlv1 
	   << G4endl;
  }

  if(eFinal <= 0.0) { 
    theParticleChange.SetMomentumChange(0.0,0.0,1.0);
    theParticleChange.SetEnergyChange(0.0);
  } else {
    theParticleChange.SetMomentumChange(nlv1.vect().unit());
    theParticleChange.SetEnergyChange(eFinal);
  }
  lv -= nlv1;
  G4double erec =  std::max(lv.e() - mass2, 0.0);
  if (verboseLevel > 1) {
    G4cout << "Recoil: " <<" m= " << mass2 << " Erec(MeV)= " << erec
	   << " 4-mom: " << lv 
	   << G4endl;
  }
 
  // the recoil is created if kinetic energy above the threshold
  if(erec > GetRecoilEnergyThreshold()) {
    G4ParticleDefinition * theDef = nullptr;
    if(Z == 1 && A == 1)       { theDef = theProton; }
    else if (Z == 1 && A == 2) { theDef = theDeuteron; }
    else if (Z == 1 && A == 3) { theDef = G4Triton::Triton(); }
    else if (Z == 2 && A == 3) { theDef = G4He3::He3(); }
    else if (Z == 2 && A == 4) { theDef = theAlpha; }
    else {
      theDef = 
	G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
    }
    G4DynamicParticle * aSec = new G4DynamicParticle(theDef, lv.vect().unit(), erec);
    theParticleChange.AddSecondary(aSec);
  } else {
    theParticleChange.SetLocalEnergyDeposit(erec);
  }

  return &theParticleChange;
}

// sample momentum transfer in the CMS system 
G4double 
G4HadronElastic::SampleInvariantT(const G4ParticleDefinition*, 
				  G4double, G4int, G4int A)
{
  static const G4double GeV2 = GeV*GeV;
  G4double tmax = pLocalTmax/GeV2;
  G4double aa, bb, cc;
  static const G4double dd = 10.;
  G4Pow* g4pow = G4Pow::GetInstance();
  if (A <= 62) {
    bb = 14.5*g4pow->Z23(A);
    aa = g4pow->powZ(A, 1.63)/bb;
    cc = 1.4*g4pow->Z13(A)/dd;
  } else {
    bb = 60.*g4pow->Z13(A);
    aa = g4pow->powZ(A, 1.33)/bb;
    cc = 0.4*g4pow->powZ(A, 0.4)/dd;
  }
  G4double q1 = 1.0 - G4Exp(-bb*tmax);
  G4double q2 = 1.0 - G4Exp(-dd*tmax);
  G4double s1 = q1*aa;
  G4double s2 = q2*cc;
  if((s1 + s2)*G4UniformRand() < s2) {
    q1 = q2;
    bb = dd;
  }
  return -GeV2*G4Log(1.0 - G4UniformRand()*q1)/bb;
}
