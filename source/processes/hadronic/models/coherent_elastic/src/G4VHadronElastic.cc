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
// $Id: G4VHadronElastic.cc,v 1.6 2010-11-19 18:50:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Geant4 Header : G4VHadronElastic
//
// Author : V.Ivanchenko 29 June 2009 (redesign old elastic model)
//  
// Modified:
//
//

#include "G4VHadronElastic.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4Pow.hh"

G4VHadronElastic::G4VHadronElastic(const G4String& name) 
  : G4HadronicInteraction(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  lowestEnergyLimit= 1.e-6*eV;  

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();
}

G4VHadronElastic::~G4VHadronElastic()
{}

G4HadFinalState* G4VHadronElastic::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();
  if(ekin <= lowestEnergyLimit) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  G4double plab = aParticle->GetTotalMomentum();

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();

  if (verboseLevel>1) {
    G4cout << "G4VHadronElastic: " 
	   << aParticle->GetDefinition()->GetParticleName() 
	   << " Plab(GeV/c)= " << plab/GeV  
	   << " Ekin(MeV) = " << ekin/MeV 
	   << " scattered off Z= " << Z 
	   << " A= " << A 
	   << G4endl;
  }
  G4ParticleDefinition * theDef = 0;

  if(Z == 1 && A == 1)       theDef = theProton;
  else if (Z == 1 && A == 2) theDef = theDeuteron;
  else if (Z == 1 && A == 3) theDef = G4Triton::Triton();
  else if (Z == 2 && A == 3) theDef = G4He3::He3();
  else if (Z == 2 && A == 4) theDef = theAlpha;
  else {
    theDef = 
      G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
  }
  G4double m2 = theDef->GetPDGMass();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,m2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  momentumCMS = p1.mag();
  G4double tmax = 4.0*momentumCMS*momentumCMS;

  // Sampling in CM system
  G4double t    = SampleInvariantT(theParticle, plab, Z, A);
  G4double phi  = G4UniformRand()*CLHEP::twopi;
  G4double cost = 1. - 2.0*t/tmax;
  G4double sint;

  // problem in sampling
  if(cost > 1.0 || cost < -1.0) {
    if(verboseLevel > 0) {
      G4cout << "G4VHadronElastic WARNING cost= " << cost
	     << " after scattering of " 
	     << aParticle->GetDefinition()->GetParticleName()
	     << " p(GeV/c)= " << plab
	     << " on " << theDef->GetParticleName()
	     << G4endl;
    }
    cost = 1.0;
    sint = 0.0;

    // normal situation
  } else  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }    
  if (verboseLevel>1) {
    G4cout << " t= " << t << " tmax= " << tmax 
	   << " Pcms= " << momentumCMS << "cos(t)=" << cost 
	   << " std::sin(t)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= momentumCMS;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),
		       std::sqrt(momentumCMS*momentumCMS + m1*m1));

  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;
  if (verboseLevel > 1) {
    G4cout <<" m= " << m1 << " Efin(MeV)= " << eFinal 
	   << " Proj: 4-mom " << lv1 << " Final: " << nlv1 
	   << G4endl;
  }
  if(eFinal <= lowestEnergyLimit) {
    if(eFinal < 0.0 && verboseLevel > 0) {
      G4cout << "G4VHadronElastic WARNING Efinal= " << eFinal
	     << " after scattering of " 
	     << aParticle->GetDefinition()->GetParticleName()
	     << " p(GeV/c)= " << plab
	     << " on " << theDef->GetParticleName()
	     << G4endl;
    }
    theParticleChange.SetEnergyChange(0.0);
    nlv1 = G4LorentzVector(0.0,0.0,0.0,m1);

  } else {
    theParticleChange.SetMomentumChange(nlv1.vect().unit());
    theParticleChange.SetEnergyChange(eFinal);
  }  

  G4LorentzVector nlv0 = lv - nlv1;
  G4double erec =  nlv0.e() - m2;
  if (verboseLevel > 1) {
    G4cout << "Recoil: " <<" m= " << m2 << " Erec(MeV)= " << erec
	   << " 4-mom: " << nlv0 
	   << G4endl;
  }
 
  if(erec > GetRecoilEnergyThreshold()) {
    G4DynamicParticle * aSec = new G4DynamicParticle(theDef, nlv0);
    theParticleChange.AddSecondary(aSec);
  } else if(erec > 0.0) {
    theParticleChange.SetLocalEnergyDeposit(erec);
  }

  return &theParticleChange;
}

// sample momentum transfer in the CMS system 
G4double 
G4VHadronElastic::SampleInvariantT(const G4ParticleDefinition* /*p*/, 
				   G4double /*ptot*/,
				   G4int /*Z*/, G4int A)
{
  static const G4double GeV2 = GeV*GeV;
  G4double tmax = 4.0*momentumCMS*momentumCMS/GeV2;
  G4double aa, bb, cc;
  G4double dd = 10.;
  G4Pow* p = G4Pow::GetInstance();
  if (A <= 62) {
    bb = 14.5*p->Z23(A);
    aa = p->powZ(A, 1.63)/bb;
    cc = 1.4*p->Z13(A)/dd;
  } else {
    bb = 60.*p->Z13(A);
    aa = p->powZ(A, 1.33)/bb;
    cc = 0.4*p->powZ(A, 0.4)/dd;
  }
  G4double q1 = 1.0 - std::exp(-bb*tmax);
  G4double q2 = 1.0 - std::exp(-dd*tmax);
  G4double s1 = q1*aa;
  G4double s2 = q2*cc;
  if((s1 + s2)*G4UniformRand() < s2) {
    q1 = q2;
    bb = dd;
  }
  return -GeV2*std::log(1.0 - G4UniformRand()*q1)/bb;
}

