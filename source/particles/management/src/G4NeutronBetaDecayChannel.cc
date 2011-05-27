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
// $Id: G4NeutronBetaDecayChannel.cc,v 1.7 2006/06/29 19:25:38 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-ref-00 $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      18 Sep  2001 H.Kurashige
//---
//      Fix energy of proton and neutrino May  2011 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4NeutronBetaDecayChannel.hh"
#include "Randomize.hh"
#include "G4RotationMatrix.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


G4NeutronBetaDecayChannel::G4NeutronBetaDecayChannel(
				       const G4String& theParentName, 
				       G4double        theBR)
                   :G4VDecayChannel("Neutron Decay"),
		    aENuCorr(-0.102)
{
  // set names for daughter particles
  if (theParentName == "neutron") {
    SetBR(theBR);
    SetParent("neutron");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "anti_nu_e");
    SetDaughter(2, "proton");
  } else if (theParentName == "anti_neutron") {
    SetBR(theBR);
    SetParent("anti_neutron");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "nu_e");
    SetDaughter(2, "anti_proton");
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4NeutronBetaDecayChannel:: constructor :";
      G4cout << " parent particle is not neutron but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

G4NeutronBetaDecayChannel::~G4NeutronBetaDecayChannel()
{
}

G4DecayProducts *G4NeutronBetaDecayChannel::DecayIt(G4double) 
{
  //  This class describes free neutron beta decay  kinemtics.
  //  This version neglects neutron/electron polarization  
  //  without Coulomb effect

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4NeutronBetaDecayChannel::DecayIt ";
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
  G4double xmax = parentmass-sumofdaughtermass;  

   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  G4double daughtermomentum[3];

  // calcurate electron energy
  G4double x;                    // Ee
  G4double p;                    // Pe
  G4double m = daughtermass[0];  //Me
  G4double w;                    // cosine of e-nu angle
  G4double r;  
  G4double r0;
  do {
      x = xmax*G4UniformRand();
      p = std::sqrt(x*(x+2.0*m));
      w = 1.0-2.0*G4UniformRand();
      r = p*(x+m)*(xmax-x)*(xmax-x)*(1.0+aENuCorr*p/(x+m)*w);
      r0 = G4UniformRand()*(xmax+m)*(xmax+m)*xmax*xmax*(1.0+aENuCorr);
  } while (r < r0);    


  //create daughter G4DynamicParticle 
  // rotation materix to lab frame
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double theta = std::acos(costheta)*rad;
  G4double phi  = twopi*G4UniformRand()*rad;
  G4RotationMatrix rm;
  rm.rotateY(theta);
  rm.rotateZ(phi);

  // daughter 0 (electron) in Z direction
  daughtermomentum[0] = p;
  G4ThreeVector direction0(0.0, 0.0, 1.0);
  direction0 = rm * direction0;
  G4DynamicParticle * daughterparticle0 
         = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle0);

  // daughter 1 (nutrino) in XZ plane
  G4double eNu;    // Enu
  eNu = (parentmass-daughtermass[2])*(parentmass+daughtermass[2])+(m*m)-2.*parentmass*(x+m);
  eNu /= 2.*(parentmass+p*w-(x+m));
  G4double cosn = w;
  G4double sinn = std::sqrt((1.0-cosn)*(1.0+cosn));

  G4ThreeVector direction1(sinn, 0.0, cosn);
  direction1 = rm * direction1;
  G4DynamicParticle * daughterparticle1 
         = new G4DynamicParticle( daughters[1], direction1*eNu);
  products->PushProducts(daughterparticle1);

  // daughter 2 (proton) at REST
  G4double eP;     // Eproton
  eP = parentmass-eNu-(x+m)-daughtermass[2];
  G4double pPx = -eNu*sinn;
  G4double pPz = -p-eNu*cosn;
  G4double pP  = std::sqrt(eP*(eP+2.*daughtermass[2]));
  G4ThreeVector direction2(pPx/pP, 0.0, pPz/pP);
    G4DynamicParticle * daughterparticle2 
         = new G4DynamicParticle( daughters[2], direction2);
  products->PushProducts(daughterparticle2);
 

 // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4NeutronBetaDecayChannel::DecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    products->DumpInfo();
  }
#endif
  return products;
}






