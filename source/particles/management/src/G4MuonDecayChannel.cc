// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuonDecayChannel.cc,v 1.4 1999-04-14 10:28:27 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      30 May  1997 H.Kurashige
//      10 June 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4MuonDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


G4MuonDecayChannel::G4MuonDecayChannel(const G4String& theParentName, 
				       G4double        theBR)
                   :G4VDecayChannel("Muon Decay",1)
{
  //#ifdef G4VERBOSE
  //if (GetVerboseLevel()>1) {
  //  G4cout << "G4MuonDecayChannel:: constructor ";
  //  G4cout << "addr[" << this << "]" << endl;
  //}
  //#endif

  // set names for daughter particles
  if (theParentName == "mu+") {
    SetBR(theBR);
    SetParent("mu+");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "nu_e");
    SetDaughter(2, "anti_nu_mu");
  } else if (theParentName == "mu-") {
    SetBR(theBR);
    SetParent("mu-");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "anti_nu_e");
    SetDaughter(2, "nu_mu");
  } else {
    //#ifdef G4VERBOSE
    //  if (GetVerboseLevel()>0) {
    //   G4cout << "G4MuonDecayChannel:: constructor :";
    //   G4cout << " parent particle is not muon but ";
    //   G4cout << theParentName << endl;
    // }
    //#endif
  }
}

G4MuonDecayChannel::~G4MuonDecayChannel()
{
}

G4DecayProducts *G4MuonDecayChannel::DecayIt(G4double) 
{
  // this version neglects muon polarization 
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for neutrinos
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4MuonDecayChannel::DecayIt ";
#endif

  if (parent == 0) FillParent();  
  if (daughters == 0) FillDaughters();
 
  // parent mass
  G4double parentmass = parent->GetPDGMass();

  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }

   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  G4double daughtermomentum[3];
  G4double momentumsum = 0.0;
  G4double energy;
  // calcurate electron energy
  G4double xmax = (1.0+daughtermass[0]*daughtermass[0]/parentmass/parentmass);
  G4double x;
  G4double r;
  do {
    do {
      r = G4UniformRand();
      x = xmax*G4UniformRand();
    } while (r < (3.0 - 2.0*x)*x*x);    
    energy = x*parentmass/2.0 - daughtermass[0];
   } while (energy <0.0);
  //create daughter G4DynamicParticle 
  // daughter 0 (electron)
  daughtermomentum[0] = sqrt(energy*energy + 2.0*energy* daughtermass[0]);
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
  G4double energy2 = parentmass*(1.0 - x/2.0); 
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
    G4cout << "G4MuonDecayChannel::DecayIt ";
    G4cout << "  create decay products in rest frame " <<endl;
    products->DumpInfo();
  }
#endif
  return products;
}






