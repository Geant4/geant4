// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DalitzDecayChannel.cc,v 1.1 1999-01-07 16:10:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

G4DalitzDecayChannel::G4DalitzDecayChannel(
			   const G4String& theParentName,
			   G4double        theBR,
			   const G4String& theLeptonName,
			   const G4String& theAntiLeptonName)
                   :G4VDecayChannel("Dalitz Decay",1)
{
  //#ifdef G4VERBOSE
  //if (GetVerboseLevel()>1) {
  //  G4cout << "G4DalitzDecayChannel:: constructor ";
  //  G4cout << "addr[" << this << "]" << endl;
  //}
  //#endif
  // set names for daughter particles
  SetParent(theParentName);
  SetBR(theBR);
  SetNumberOfDaughters(3);
  SetDaughter(idGamma, "gamma");
  SetDaughter(idLepton, theLeptonName);
  SetDaughter(idAntiLepton, theAntiLeptonName);
  //
  //#ifdef G4VERBOSE
  //if (GetVerboseLevel()>1) DumpInfo();
  //#endif
}

G4DalitzDecayChannel::~G4DalitzDecayChannel()
{
}

G4DecayProducts *G4DalitzDecayChannel::DecayIt(G4double) 
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4DalitzDecayChannel::DecayIt ";
#endif 
  if (parent == NULL) FillParent();  
  if (daughters == NULL) FillDaughters();

  // parent mass
  G4double parentmass = parent->GetPDGMass();
 
 //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
 
  //daughters'mass
  G4double leptonmass = daughters[idLepton]->GetPDGMass(); 

 // Generate t ( = exp(x):mass Square of (l+ l-) system) 
  G4double xmin  = 2.0*log(2.0*leptonmass);
  G4double xmax  = 2.0*log(parentmass);
  G4double wmax = 1.5;
  G4double x, w, ww, w1, w2, w3, t;
  do {
    x = G4UniformRand()*(xmax-xmin) + xmin;
    w = G4UniformRand()*wmax;
    t = exp(x);
    w1 = (1.0-4.0*leptonmass*leptonmass/t);
    if ( w1 > 0.0) {
      w2 = ( 1.0 + 2.0*leptonmass*leptonmass/t );
      w3 = ( 1.0 - t/parentmass/parentmass );
      w3 = w3 * w3 * w3;
      ww = w3 * w2 * sqrt(w1);
    } else {
      ww = 0.0;
    }
  } while (w > ww);    
 
  // calculate gamma momentum
  G4double Pgamma = 
      G4PhaseSpaceDecayChannel::Pmx(parentmass, 0.0, sqrt(t)); 
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi  = 2.0*M_PI*G4UniformRand()*rad;
  G4ParticleMomentum gdirection(sintheta*cos(phi),sintheta*sin(phi),costheta);

  //create G4DynamicParticle for gamma 
  G4DynamicParticle * gammaparticle
      = new G4DynamicParticle(daughters[idGamma] , gdirection, Pgamma);

  // calcurate beta of (l+ l-)system
  G4double beta = Pgamma/(parentmass-Pgamma);

  // calculate lepton momentum in the rest frame of (l+ l-)system
  G4double Plepton = 
      G4PhaseSpaceDecayChannel::Pmx(sqrt(t),leptonmass, leptonmass); 
  G4double Elepton = sqrt(Plepton*Plepton + leptonmass*leptonmass );
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = sqrt((1.0 - costheta)*(1.0 + costheta));
  phi  = 2.0*M_PI*G4UniformRand()*rad;
  G4ParticleMomentum ldirection(sintheta*cos(phi),sintheta*sin(phi),costheta);
  //create G4DynamicParticle for leptons  in the rest frame of (l+ l-)system
  G4DynamicParticle * leptonparticle 
    = new G4DynamicParticle(daughters[idLepton] , 
			    ldirection, Elepton-leptonmass );
  G4DynamicParticle * antileptonparticle 
    = new G4DynamicParticle(daughters[idAntiLepton] , 
			    -1.0*ldirection, Elepton-leptonmass );
  //boost leptons in the rest frame of the parent 
  G4LorentzVector p4 = leptonparticle->Get4Momentum();
  p4.boost( -1.0*gdirection.x()*beta, -1.0*gdirection.y()*beta, -1.0*gdirection.z()*beta);
  leptonparticle->Set4Momentum(p4);
  p4 = antileptonparticle->Get4Momentum();
  p4.boost( -1.0*gdirection.x()*beta, -1.0*gdirection.y()*beta, -1.0*gdirection.z()*beta);
  antileptonparticle->Set4Momentum(p4);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;
  products->PushProducts(gammaparticle);
  products->PushProducts(leptonparticle);
  products->PushProducts(antileptonparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4DalitzDecayChannel::DecayIt ";
     G4cout << "  create decay products in rest frame " <<endl;
     products->DumpInfo();
  }
#endif
  return products;
}





