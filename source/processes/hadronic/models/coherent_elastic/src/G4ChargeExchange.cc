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
// $Id: G4ChargeExchange.cc 91897 2015-08-10 09:55:12Z gcosmo $
//
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


G4ChargeExchange::G4ChargeExchange() : G4HadronicInteraction("Charge Exchange")
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );

  lowEnergyRecoilLimit = 100.*keV;  
  lowestEnergyLimit    = 1.*MeV;  

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
  theHe3      = G4He3::He3();
}

G4ChargeExchange::~G4ChargeExchange()
{}

G4HadFinalState* G4ChargeExchange::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  if(ekin <= lowestEnergyLimit || A < 3) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4double plab = aParticle->GetTotalMomentum();

  if (verboseLevel > 1)
    G4cout << "G4ChargeExchange::DoIt: Incident particle plab="
	   << plab/GeV << " GeV/c "
	   << " ekin(MeV) = " << ekin/MeV << "  "
	   << aParticle->GetDefinition()->GetParticleName() << G4endl;

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4int N = A - Z;
  G4int projPDG = theParticle->GetPDGEncoding();
  if (verboseLevel > 1)
    G4cout << "G4ChargeExchange for " << theParticle->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z
	   << " A= " << A << " N= " << N
	   << G4endl;

  G4ParticleDefinition * theDef = 0;

  G4double mass2 = G4NucleiProperties::GetNuclearMass((G4double)A, (G4double)Z);
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv0(0.0,0.0,0.0,mass2);

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
    if(G4UniformRand()*A < G4double(Z)) {
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
    if(G4UniformRand()*A < G4double(Z)) {
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
  else {
    theDef = 
      G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
  }
  if(!theSecondary) { return &theParticleChange; }

  G4double m11 = theSecondary->GetPDGMass();
  G4double m21 = theDef->GetPDGMass();
  if(theRecoil)  { m21 += theRecoil->GetPDGMass(); }
  else           { theRecoil = theDef; }

  G4double etot = lv0.e() + lv1.e();

  // kinematiacally impossible
  if(etot < m11 + m21) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4ThreeVector p1 = lv1.vect();
  G4double e1 = 0.5*etot*(1.0 - (m21*m21 - m11*m11)/(etot*etot));
  //  G4double e2 = etot - e1;
  G4double ptot = std::sqrt(e1*e1 - m11*m11);

  G4double tmax = 4.0*ptot*ptot;
  G4double g2 = GeV*GeV; 

  G4double t = g2*SampleT(tmax/g2, A);

  if(verboseLevel>1)
    G4cout <<"## G4ChargeExchange t= " << t << " tmax= " << tmax
	   << " ptot= " << ptot << G4endl;

  // Sampling in CM system
  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  if(std::abs(cost) > 1.0) cost = 1.0;
  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  //if (verboseLevel > 1)
  //  G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;

  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),e1);
  G4LorentzVector nlv0 = lv0 + lv1 - nlv1;

  nlv0.boost(bst);
  nlv1.boost(bst);

  theParticleChange.SetStatusChange(stopAndKill);
  theParticleChange.SetEnergyChange(0.0);
  G4DynamicParticle * aSec = new G4DynamicParticle(theSecondary, nlv1);
  theParticleChange.AddSecondary(aSec);

  G4double erec =  nlv0.e() - m21;

  //G4cout << "erec= " <<erec << " Esec= " << aSec->GetKineticEnergy() << G4endl;  

  if(theHyperon) {
    theParticleChange.SetLocalEnergyDeposit(erec);
    aSec = new G4DynamicParticle();
    aSec->SetDefinition(theRecoil);
    aSec->SetKineticEnergy(0.0);
  } else if(erec > lowEnergyRecoilLimit) {
    aSec = new G4DynamicParticle(theRecoil, nlv0);
    theParticleChange.AddSecondary(aSec);
  } else {
    if(erec < 0.0) erec = 0.0;
    theParticleChange.SetLocalEnergyDeposit(erec);
  }
  return &theParticleChange;
}

G4double G4ChargeExchange::SampleT(G4double tmax, G4double A)
{
  G4double aa, bb, cc, dd;
  if (A <= 62.) {
    aa = G4Pow::GetInstance()->powA(A, 1.63);
    bb = 14.5*G4Pow::GetInstance()->powA(A, 0.66);
    cc = 1.4*G4Pow::GetInstance()->powA(A, 0.33);
    dd = 10.;
  } else {
    aa = G4Pow::GetInstance()->powA(A, 1.33);
    bb = 60.*G4Pow::GetInstance()->powA(A, 0.33);
    cc = 0.4*G4Pow::GetInstance()->powA(A, 0.40);
    dd = 10.;
  }
  G4double x1 = (1.0 - G4Exp(-tmax*bb))*aa/bb;
  G4double x2 = (1.0 - G4Exp(-tmax*dd))*cc/dd;
  
  G4double t;
  G4double y = bb;
  if(G4UniformRand()*(x1 + x2) < x2) y = dd;

  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = 0;
  do {
    t = -G4Log(G4UniformRand())/y;
  } while ( (t > tmax) &&
            ++loopCounter < maxNumberOfLoops );  /* Loop checking, 10.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    t = 0.0;
  }

  return t;
}

