//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TauLeptonicDecayChannel.cc,v 1.1 2002-03-08 08:47:55 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May  1997 H.Kurashige
//
//      Fix bug in calcuration of electron energy in DecayIt 28 Feb. 01 H.Kurashige 
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4TauLeptonicDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


G4TauLeptonicDecayChannel::G4TauLeptonicDecayChannel(
						     const G4String& theParentName, 
						     G4double        theBR,
						     const G4String& theLeptonName)
                   :G4VDecayChannel("Tau Leptonic Decay",1)
{
  // set names for daughter particles
  if (theParentName == "tau+") {
    SetBR(theBR);
    SetParent("tau+");
    SetNumberOfDaughters(3);
    if (theLeptonName[0]=='e'){
      SetDaughter(0, "e+");
      SetDaughter(1, "nu_e");
      SetDaughter(2, "anti_nu_tau");
    } else {
      SetDaughter(0, "mu+");
      SetDaughter(1, "nu_mu");
      SetDaughter(2, "anti_nu_tau");
    } 
  } else if (theParentName == "tau-") {
    SetBR(theBR);
    SetParent("tau-");
    SetNumberOfDaughters(3);
    if (theLeptonName[0]=='e'){
      SetDaughter(0, "e-");
      SetDaughter(1, "anti_nu_e");
      SetDaughter(2, "nu_tau");
    } else {
      SetDaughter(0, "mu-");
      SetDaughter(1, "anti_nu_mu");
      SetDaughter(2, "nu_tau");
    } 
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4TauLeptonicDecayChannel:: constructor :";
      G4cout << " parent particle is not tau but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

G4TauLeptonicDecayChannel::~G4TauLeptonicDecayChannel()
{
}

G4DecayProducts *G4TauLeptonicDecayChannel::DecayIt(G4double) 
{
  // this version neglects muon polarization 
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for neutrinos
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4TauLeptonicDecayChannel::DecayIt ";
#endif

  if (parent == 0) FillParent();  
  if (daughters == 0) FillDaughters();
 
  // parent mass
  G4double parentmass = parent->GetPDGMass();

  //daughters'mass
  G4double daughtermass[3]; 
  for (G4int index=0; index<3; index++){
    daughtermass[index] = daughters[index]->GetPDGMass();
  }

   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  G4double daughtermomentum[3];
  G4double energy;

  // calcurate lepton momentum
  G4double pmax = (parentmass*parentmass-daughtermass[0]*daughtermass[0])/2./parentmass;
  G4double p, e;
  G4double r;
  do {
    // determine momentum/energy
    r = G4UniformRand();
    p = pmax*G4UniformRand();
    e = sqrt(p*p + daughtermass[0]*daughtermass[0]);
  } while (r > spectrum(p,e,parentmass,daughtermass[0]) );
  energy = e- daughtermass[0];

  //create daughter G4DynamicParticle 
  // daughter 0 (lepton)
  daughtermomentum[0] = p;
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = sqrt((1.0-costheta)*(1.0+costheta));
  phi  = 2.0*M_PI*G4UniformRand()*rad;
  sinphi = sin(phi);
  cosphi = cos(phi);
  G4ThreeVector direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle 
         = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  // daughter 1 ,2 (nutrinos)
  // create neutrinos in the C.M frame of two neutrinos
  G4double energy2 = parentmass-e; 
  G4double vmass   = sqrt((energy2-daughtermomentum[0])*(energy2+daughtermomentum[0]));
  G4double beta = -1.0*daughtermomentum[0]/energy2;
  G4double costhetan = 2.*G4UniformRand()-1.0;
  G4double sinthetan = sqrt((1.0-costhetan)*(1.0+costhetan));
  G4double phin  = 2.0*M_PI*G4UniformRand()*rad;
  G4double sinphin = sin(phin);
  G4double cosphin = cos(phin);

  G4ThreeVector direction1(sinthetan*cosphin,sinthetan*sinphin,costhetan);
  G4DynamicParticle * daughterparticle1 
         = new G4DynamicParticle( daughters[1], direction1*(vmass/2.));
  G4DynamicParticle * daughterparticle2
         = new G4DynamicParticle( daughters[2], direction1*(-1.0*vmass/2.));

  // boost to the muon rest frame
  G4LorentzVector p4;
  p4 = daughterparticle1->Get4Momentum();
  p4.boost( direction0.x()*beta, direction0.y()*beta, direction0.z()*beta);
  daughterparticle1->Set4Momentum(p4);
  p4 = daughterparticle2->Get4Momentum();
  p4.boost( direction0.x()*beta, direction0.y()*beta, direction0.z()*beta);
  daughterparticle2->Set4Momentum(p4);
  products->PushProducts(daughterparticle1);
  products->PushProducts(daughterparticle2);
  daughtermomentum[1] = daughterparticle1->GetTotalMomentum();
  daughtermomentum[2] = daughterparticle2->GetTotalMomentum();


 // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4TauLeptonicDecayChannel::DecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    products->DumpInfo();
  }
#endif
  return products;
}




G4double G4TauLeptonicDecayChannel::spectrum(G4double p,
					     G4double e,
					     G4double mtau,
					     G4double ml)
{
  G4double f1;
  f1 = 3.0*e*(mtau*mtau+ml*ml)-4.0*mtau*e*e-2.0*mtau*ml*ml;
  return p*(f1)/(mtau*mtau*mtau*mtau)/(0.6);
}



