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
// $Id: G4DalitzDecayChannel.cc 95906 2016-03-02 10:56:50Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

G4DalitzDecayChannel::G4DalitzDecayChannel()
  :G4VDecayChannel()
{
}

G4DalitzDecayChannel::G4DalitzDecayChannel(
			   const G4String& theParentName,
			   G4double        theBR,
			   const G4String& theLeptonName,
			   const G4String& theAntiLeptonName)
                   :G4VDecayChannel("Dalitz Decay",1)
{
  // set names for daughter particles
  SetParent(theParentName);
  SetBR(theBR);
  SetNumberOfDaughters(3);
  G4String gammaName = "gamma";
  SetDaughter(idGamma, gammaName);
  SetDaughter(idLepton, theLeptonName);
  SetDaughter(idAntiLepton, theAntiLeptonName);
}

G4DalitzDecayChannel::~G4DalitzDecayChannel()
{
}

G4DalitzDecayChannel::G4DalitzDecayChannel(const G4DalitzDecayChannel &right):
  G4VDecayChannel(right)
{
}

G4DalitzDecayChannel & G4DalitzDecayChannel::operator=(const G4DalitzDecayChannel & right)
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

G4DecayProducts *G4DalitzDecayChannel::DecayIt(G4double) 
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4DalitzDecayChannel::DecayIt ";
#endif 
  CheckAndFillParent();
  CheckAndFillDaughters();

  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();
 
 //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);
 
  //daughters'mass
  G4double leptonmass = G4MT_daughters[idLepton]->GetPDGMass(); 

 // Generate t ( = std::exp(x):mass Square of (l+ l-) system) 
  G4double xmin  = 2.0*std::log(2.0*leptonmass);
  G4double xmax  = 2.0*std::log(parentmass);
  G4double wmax = 1.5;
  G4double x, w, ww, w1, w2, w3, t;
  const size_t MAX_LOOP = 10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    x = G4UniformRand()*(xmax-xmin) + xmin;
    w = G4UniformRand()*wmax;
    t = std::exp(x);
    w1 = (1.0-4.0*leptonmass*leptonmass/t);
    if ( w1 > 0.0) {
      w2 = ( 1.0 + 2.0*leptonmass*leptonmass/t );
      w3 = ( 1.0 - t/parentmass/parentmass );
      w3 = w3 * w3 * w3;
      ww = w3 * w2 * std::sqrt(w1);
    } else {
      ww = 0.0;
    }
    if (w <= ww) break;    
  }

  // calculate gamma momentum
  G4double Pgamma = 
      G4PhaseSpaceDecayChannel::Pmx(parentmass, 0.0, std::sqrt(t)); 
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector gdirection(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);

  //create G4DynamicParticle for gamma 
  G4DynamicParticle * gammaparticle
      = new G4DynamicParticle(G4MT_daughters[idGamma] , gdirection, Pgamma);

  // calcurate beta of (l+ l-)system
  G4double beta = Pgamma/(parentmass-Pgamma);

  // calculate lepton momentum in the rest frame of (l+ l-)system
  G4double Plepton = 
      G4PhaseSpaceDecayChannel::Pmx(std::sqrt(t),leptonmass, leptonmass); 
  G4double Elepton = std::sqrt(Plepton*Plepton + leptonmass*leptonmass );
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector ldirection(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);
  //create G4DynamicParticle for leptons  in the rest frame of (l+ l-)system
  G4DynamicParticle * leptonparticle 
    = new G4DynamicParticle(G4MT_daughters[idLepton] , 
			    ldirection, Elepton-leptonmass );
  G4DynamicParticle * antileptonparticle 
    = new G4DynamicParticle(G4MT_daughters[idAntiLepton] , 
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
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}





