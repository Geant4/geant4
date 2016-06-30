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
// $Id: G4TauLeptonicDecayChannel.cc 95906 2016-03-02 10:56:50Z gcosmo $
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4TauLeptonicDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


G4TauLeptonicDecayChannel::G4TauLeptonicDecayChannel()
  :G4VDecayChannel()
{
}


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
    if ((theLeptonName=="e-"||theLeptonName=="e+")){
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
    if ((theLeptonName=="e-"||theLeptonName=="e+")){
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

G4TauLeptonicDecayChannel::G4TauLeptonicDecayChannel(const G4TauLeptonicDecayChannel &right):
  G4VDecayChannel(right)
{
}

G4TauLeptonicDecayChannel & G4TauLeptonicDecayChannel::operator=(const G4TauLeptonicDecayChannel & right)
{
  if (this != &right) { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters >0 ) {
      if (daughters_name !=0) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      //copy daughters name
      for (G4int index=0; index < numberOfDaughters; index++) {
          daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }
  return *this;
}

G4DecayProducts *G4TauLeptonicDecayChannel::DecayIt(G4double) 
{
  // this version neglects muon polarization 
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for neutrinos
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4TauLeptonicDecayChannel::DecayIt ";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();
 
  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();

  //daughters'mass
  const G4int N_DAUGHTER=3;
  G4double daughtermass[N_DAUGHTER]; 
  for (G4int index=0; index<N_DAUGHTER; index++){
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
  }

   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  G4double daughtermomentum[N_DAUGHTER];

  // calcurate lepton momentum
  G4double pmax = (parentmass*parentmass-daughtermass[0]*daughtermass[0])/2./parentmass;
  G4double p, e;
  G4double r;
  const size_t MAX_LOOP=10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    // determine momentum/energy
    r = G4UniformRand();
    p = pmax*G4UniformRand();
    e = std::sqrt(p*p + daughtermass[0]*daughtermass[0]);
    if (r < spectrum(p,e,parentmass,daughtermass[0]) ) break;
  } 

  //create daughter G4DynamicParticle 
  // daughter 0 (lepton)
  daughtermomentum[0] = p;
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
  G4ThreeVector direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle 
         = new G4DynamicParticle( G4MT_daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  // daughter 1 ,2 (nutrinos)
  // create neutrinos in the C.M frame of two neutrinos
  G4double energy2 = parentmass-e; 
  G4double vmass   = std::sqrt((energy2-daughtermomentum[0])*(energy2+daughtermomentum[0]));
  G4double beta = -1.0*daughtermomentum[0]/energy2;
  G4double costhetan = 2.*G4UniformRand()-1.0;
  G4double sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
  G4double phin  = twopi*G4UniformRand()*rad;
  G4double sinphin = std::sin(phin);
  G4double cosphin = std::cos(phin);

  G4ThreeVector direction1(sinthetan*cosphin,sinthetan*sinphin,costhetan);
  G4DynamicParticle * daughterparticle1 
         = new G4DynamicParticle( G4MT_daughters[1], direction1*(vmass/2.));
  G4DynamicParticle * daughterparticle2
         = new G4DynamicParticle( G4MT_daughters[2], direction1*(-1.0*vmass/2.));

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



