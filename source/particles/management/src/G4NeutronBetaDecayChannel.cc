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
// $Id: G4NeutronBetaDecayChannel.cc 95906 2016-03-02 10:56:50Z gcosmo $
// GEANT4 tag $Name: geant4-09-04-ref-00 $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      18 Sep  2001 H.Kurashige
//---
//      Fix energy of proton and neutrino     May 2011 H.Kurashige
//      Fix direction of proton and neutrino  Nov 2013 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4NeutronBetaDecayChannel.hh"
#include "Randomize.hh"
#include "G4RotationMatrix.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

G4NeutronBetaDecayChannel::G4NeutronBetaDecayChannel()
                   :G4VDecayChannel(),
		    aENuCorr(-0.102)
{
}

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

G4NeutronBetaDecayChannel::G4NeutronBetaDecayChannel(const G4NeutronBetaDecayChannel &right)
    : G4VDecayChannel(right),
      aENuCorr(-0.102)
{
}


G4NeutronBetaDecayChannel & G4NeutronBetaDecayChannel::operator=(const G4NeutronBetaDecayChannel & right)
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

G4DecayProducts *G4NeutronBetaDecayChannel::DecayIt(G4double) 
{
  //  This class describes free neutron beta decay  kinemtics.
  //  This version neglects neutron/electron polarization  
  //  without Coulomb effect

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4NeutronBetaDecayChannel::DecayIt ";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();
 
  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();

  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }
  G4double xmax = parentmass-sumofdaughtermass;  

   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  G4double daughtermomentum[3];

  // calcurate electron energy
  G4double x;                    // Ee
  G4double p;                    // Pe
  G4double dm = daughtermass[0]; //Me
  G4double w;                    // cosine of e-nu angle
  G4double r;  
  G4double r0;
  const size_t MAX_LOOP=10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
      x = xmax*G4UniformRand();
      p = std::sqrt(x*(x+2.0*dm));
      w = 1.0-2.0*G4UniformRand();
      r = p*(x+dm)*(xmax-x)*(xmax-x)*(1.0+aENuCorr*p/(x+dm)*w);
      r0 = G4UniformRand()*(xmax+dm)*(xmax+dm)*xmax*xmax*(1.0+aENuCorr);
      if  (r > r0) break;
   }  

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
         = new G4DynamicParticle( G4MT_daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle0);

  // daughter 1 (nutrino) in XZ plane
  G4double eNu;    // Enu
  eNu = (parentmass-daughtermass[2])*(parentmass+daughtermass[2])+(dm*dm)-2.*parentmass*(x+dm);
  eNu /= 2.*(parentmass+p*w-(x+dm));
  G4double cosn = w;
  G4double phin  = twopi*G4UniformRand()*rad;
  G4double sinn = std::sqrt((1.0-cosn)*(1.0+cosn));

  G4ThreeVector direction1(sinn*std::cos(phin), sinn*std::sin(phin), cosn);
  direction1 = rm * direction1;
  G4DynamicParticle * daughterparticle1 
         = new G4DynamicParticle( G4MT_daughters[1], direction1*eNu);
  products->PushProducts(daughterparticle1);

  // daughter 2 (proton) at REST
  G4double eP;     // Eproton
  eP = parentmass-eNu-(x+dm)-daughtermass[2];
  G4double pPx = -eNu*sinn;
  G4double pPz = -p-eNu*cosn;
  G4double pP  = std::sqrt(eP*(eP+2.*daughtermass[2]));
  G4ThreeVector direction2(pPx/pP*std::cos(phin), pPx/pP*std::sin(phin), pPz/pP);
    G4DynamicParticle * daughterparticle2 
         = new G4DynamicParticle( G4MT_daughters[2], direction2);
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






