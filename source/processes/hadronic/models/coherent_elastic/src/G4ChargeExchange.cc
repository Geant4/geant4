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
// $Id: G4ChargeExchange.cc,v 1.3 2006/06/29 20:09:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//
// G4 Model: Charge and strangness exchange based on G4LightMedia model
//           28 May 2006 V.Ivanchenko
//
// Modified:
//

#include "G4ChargeExchange.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4QElasticCrossSection.hh"
#include "G4VQCrossSection.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "Randomize.hh"
#include "G4HadronElastic.hh"


G4ChargeExchange::G4ChargeExchange(G4HadronElastic* hel, G4double elim,
                                   G4double plow, G4double ehigh)
: G4HadronicInteraction(),
  fElastic(hel),
  native(false),
  ekinlim(elim),
  plablow(plow),
  ekinhigh(ehigh)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( DBL_MAX );
  verboseLevel= 0;
  if(!fElastic) {
    native = true;
    fElastic = new G4HadronElastic(elim, plow, ehigh);
  }
  qCManager   = fElastic->GetCS();
  hElastic    = fElastic->GetHElastic();

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theAProton  = G4AntiProton::AntiProton();
  theANeutron = G4AntiNeutron::AntiNeutron();
  thePiPlus   = G4PionPlus::PionPlus();
  thePiMinus  = G4PionMinus::PionMinus();
  thePiZero   = G4PionZero::PionZero();
  theKPlus    = G4KaonPlus::KaonPlus();
  theKMinus   = G4KaonMinus::KaonMinus();
  theK0S      = G4KaonZeroShort::KaonZeroShort();
  theK0L      = G4KaonZeroLong::KaonZeroLong();
  theL        = G4Lambda::Lambda();
  theAntiL    = G4AntiLambda::AntiLambda();
  theSPlus    = G4SigmaPlus::SigmaPlus();
  theASPlus   = G4AntiSigmaPlus::AntiSigmaPlus();
  theSMinus   = G4SigmaMinus::SigmaMinus();
  theASMinus  = G4AntiSigmaMinus::AntiSigmaMinus();
  theS0       = G4SigmaZero::SigmaZero();
  theAS0      = G4AntiSigmaZero::AntiSigmaZero();
  theXiMinus  = G4XiMinus::XiMinus();
  theXi0      = G4XiZero::XiZero();
  theAXiMinus = G4AntiXiMinus::AntiXiMinus();
  theAXi0     = G4AntiXiZero::AntiXiZero();
  theOmega    = G4OmegaMinus::OmegaMinus();
  theAOmega   = G4AntiOmegaMinus::AntiOmegaMinus();
  theD        = G4Deuteron::Deuteron();
  theT        = G4Triton::Triton();
  theA        = G4Alpha::Alpha();
  theA        = G4He3::He3();
}

G4ChargeExchange::~G4ChargeExchange()
{
  if(native) delete fElastic;
}

G4HadFinalState* G4ChargeExchange::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  const G4HadProjectile* aParticle = &aTrack;
  G4double aTarget = targetNucleus.GetN();
  G4double zTarget = targetNucleus.GetZ();
  theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
  theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
  G4int Z = static_cast<G4int>(zTarget);
  G4int A = static_cast<G4int>(aTarget);
  if(A < 3) return &theParticleChange;

  G4double plab = aParticle->GetTotalMomentum();
  G4double ekin = aParticle->GetKineticEnergy();
  if (verboseLevel > 1)
    G4cout << "G4ChargeExchange::DoIt: Incident particle plab="
	   << plab/GeV << " GeV/c "
	   << " ekin(MeV) = " << ekin/MeV << "  "
	   << aParticle->GetDefinition()->GetParticleName() << G4endl;

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();

  G4int N = A - Z;
  G4int projPDG = theParticle->GetPDGEncoding();
  if (verboseLevel > 1)
    G4cout << "G4ChargeExchange for " << theParticle->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z
	   << " A= " << A << " N= " << N
	   << G4endl;

  G4ParticleDefinition * theDef = 0;

  if (Z == 1 && A == 3) theDef = theT;
  else if (Z == 2 && A == 3) theDef = theHe3;
  else if (Z == 2 && A == 4) theDef = theA;
  else theDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);

  G4double m2 = theDef->GetPDGMass();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv0(0.0,0.0,0.0,m2);

  G4LorentzVector lv  = lv0 + lv1;
  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);
  lv0.boost(-bst);

  // Sample final particles
  G4bool theHyperon = false;
  G4ParticleDefinition* theRecoil = 0;
  G4ParticleDefinition* theSecondary = 0;

  if(theParticle == theProton) {
    theSecondary = theNeutron;
    Z++;
  } else if(theParticle == theNeutron) {
    theSecondary = theProton;
    Z--;
  } else if(theParticle == thePiPlus) {
    theSecondary = thePiZero;
    Z++;
  } else if(theParticle == thePiMinus) {
    theSecondary = thePiZero;
    Z--;
  } else if(theParticle == theKPlus) {
    if(G4UniformRand()<0.5) theSecondary = theK0S;
    else  theSecondary = theK0L;
    Z++;
  } else if(theParticle == theKMinus) {
    if(G4UniformRand()<0.5) theSecondary = theK0S;
    else  theSecondary = theK0L;
    Z--;
  } else if(theParticle == theK0S || theParticle == theK0L) {
    if(G4UniformRand()*aTarget < zTarget) {
      theSecondary = theKPlus;
      Z--;
    } else {
      theSecondary = theKMinus;
      Z++;
    }
  } else if(theParticle == theANeutron) {
    theSecondary = theAProton;
    Z++;
  } else if(theParticle == theAProton) {
    theSecondary = theANeutron;
    Z--;
  } else if(theParticle == theL) {
    G4double x = G4UniformRand();
    if(G4UniformRand()*aTarget < zTarget) {
      if(x < 0.2) {
        theSecondary = theS0;
      } else if (x < 0.4) {
        theSecondary = theSPlus;
        Z--;
      } else if (x < 0.6) {
        theSecondary = theProton;
	theRecoil = theL;
        theHyperon = true;
	A--;
      } else if (x < 0.8) {
        theSecondary = theProton;
	theRecoil = theS0;
        theHyperon = true;
        A--;
      } else {
        theSecondary = theNeutron;
	theRecoil = theSPlus;
        theHyperon = true;
        A--;
      }
    } else {
      if(x < 0.2) {
        theSecondary = theS0;
      } else if (x < 0.4) {
        theSecondary = theSMinus;
        Z++;
      } else if (x < 0.6) {
        theSecondary = theNeutron;
	theRecoil = theL;
        A--;
        theHyperon = true;
      } else if (x < 0.8) {
        theSecondary = theNeutron;
	theRecoil = theS0;
        theHyperon = true;
        A--;
      } else {
        theSecondary = theProton;
	theRecoil = theSMinus;
        theHyperon = true;
        A--;
      }
    }
  }

  if (Z == 1 && A == 2) theDef = theD;
  else if (Z == 1 && A == 3) theDef = theT;
  else if (Z == 2 && A == 3) theDef = theHe3;
  else if (Z == 2 && A == 4) theDef = theA;
  else theDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);

  G4double m11 = theSecondary->GetPDGMass();
  G4double m21 = theDef->GetPDGMass();
  if(theRecoil)  m21 += theRecoil->GetPDGMass();
  else           theRecoil = theDef;

  G4double etot = lv0.e() + lv1.e();
  if(etot < m11 + m21) return &theParticleChange;

  G4ThreeVector p1 = lv1.vect();
  G4double e1 = 0.5*etot*(1.0 + (m21*m21 - m11*m11)/(etot*etot));
  G4double e2 = etot - e1;
  G4double ptot = std::sqrt(e1*e1 - m11*m11);

  G4double tmax = 4.0*ptot*ptot;
  G4double t = 0.0;

  // Choose generator
  G4ElasticGenerator gtype = fLElastic;
  if ((theParticle == theProton || theParticle == theNeutron) && Z <= 2) {
    gtype = fQElastic;
    if(Z == 1 && N == 2) N = 1;
    else if (Z == 2 && N == 1) N = 2;
  } else if(ekin >= ekinhigh) {
    gtype = fHElastic;
  } else if(plab <= plablow) {
    gtype = fSWave;
  }

  // Sample t
  if(gtype == fQElastic) {
    if (verboseLevel > 1)
      G4cout << "G4ChargeExchange: Z= " << Z << " N= "
	     << N << " pdg= " <<  projPDG
	     << " mom(GeV)= " << plab/GeV << "  " << qCManager << G4endl;
    G4double cs = qCManager->GetCrossSection(false,plab,Z,N,projPDG);
    if(cs > 0.0) t = qCManager->GetExchangeT(Z,N,projPDG);
    else gtype = fSWave;
  }

  if(gtype == fSWave)         t = G4UniformRand()*tmax;
  else if(gtype == fHElastic) t = hElastic->SampleT(theParticle,plab,Z,A);
  else if(gtype == fLElastic) t = GeV*GeV*fElastic->SampleT(ptot,m1,m2,aTarget);

  if(verboseLevel>1)
    G4cout <<"type= " << gtype <<" t= " << t << " tmax= " << tmax
	   << " ptot= " << ptot << G4endl;

  // Sampling in CM system
  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  if(std::abs(cost) > 1.0) cost = -1.0 + 2.0*G4UniformRand();
  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  if (verboseLevel > 1)
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;

  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  p1 = p1.unit();
  v1.rotateUz(p1);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),e1);
  G4LorentzVector nlv0(-v1.x(),-v1.y(),-v1.z(),e2);

  nlv0.boost(bst);
  nlv1.boost(bst);

  theParticleChange.SetStatusChange(stopAndKill);
  G4DynamicParticle * aSec = new G4DynamicParticle(theSecondary, nlv1);
  theParticleChange.AddSecondary(aSec);

  G4double erec =  nlv0.e() - m21;
  if(theHyperon) {
    theParticleChange.SetLocalEnergyDeposit(erec);
    aSec = new G4DynamicParticle();
    aSec->SetDefinition(theRecoil);
    aSec->SetKineticEnergy(0.0);
  } else if(erec > ekinlim) {
    aSec = new G4DynamicParticle(theRecoil, nlv0);
    theParticleChange.AddSecondary(aSec);
  } else {
    theParticleChange.SetLocalEnergyDeposit(erec);
  }
  return &theParticleChange;
}
