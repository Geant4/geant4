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
// G4MuonDecayChannel class implementation
//
// Author: H.Kurashige, 30 May 1997
// Contributions:
// - 2005 - M.Melissas, J.Brunner CPPM/IN2P3
//   Added V-A fluxes for neutrinos using a new algorithm, 2005
// --------------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4MuonDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"

G4MuonDecayChannel::G4MuonDecayChannel()
  : G4VDecayChannel()
{
}

G4MuonDecayChannel::G4MuonDecayChannel(const G4String& theParentName, 
                                             G4double theBR)
  : G4VDecayChannel("Muon Decay", 1)
{
  // set names for daughter particles
  if (theParentName == "mu+")
  {
    SetBR(theBR);
    SetParent("mu+");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "nu_e");
    SetDaughter(2, "anti_nu_mu");
  }
  else if (theParentName == "mu-")
  {
    SetBR(theBR);
    SetParent("mu-");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "anti_nu_e");
    SetDaughter(2, "nu_mu");
  }
  else
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4MuonDecayChannel:: constructor :";
      G4cout << " parent particle is not muon but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

G4MuonDecayChannel::G4MuonDecayChannel(const G4MuonDecayChannel& right)
  : G4VDecayChannel(right)
{
}

G4MuonDecayChannel::~G4MuonDecayChannel()
{
}

G4MuonDecayChannel&
G4MuonDecayChannel::operator=(const G4MuonDecayChannel& right)
{
  if (this != &right)
  { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters >0 )
    {
      if (daughters_name !=0) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      // copy daughters name
      for (G4int index=0; index < numberOfDaughters; ++index)
      {
        daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }
  return *this;
}

G4DecayProducts* G4MuonDecayChannel::DecayIt(G4double) 
{
  // this version neglects muon polarization,and electron mass  
  //              assumes the pure V-A coupling
  //              the Neutrinos are correctly V-A

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4MuonDecayChannel::DecayIt ";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();
 
  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();
  const G4int N_DAUGHTER=3;

  // daughters'mass
  G4double daughtermass[N_DAUGHTER]; 
  //G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<N_DAUGHTER; ++index)
  {
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    //sumofdaughtermass += daughtermass[index];
  }

  // create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle* parentparticle
    = new G4DynamicParticle( G4MT_parent, dummy, 0.0);
  // create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  G4double daughtermomentum[N_DAUGHTER];
  // calculate electron energy
  G4double xmax = (1.0+daughtermass[0]*daughtermass[0]/parentmass/parentmass);
  G4double x;
  
  G4double Ee,Ene;
  
  G4double gam;
  G4double EMax=parentmass/2-daughtermass[0];
  
  const std::size_t MAX_LOOP=1000;
   //Generating Random Energy
  for (std::size_t loop1=0; loop1<MAX_LOOP; ++loop1)
  {
    Ee=G4UniformRand();
    for (std::size_t loop2 =0; loop2<MAX_LOOP; ++loop2)
    {
      x=xmax*G4UniformRand();
      gam=G4UniformRand();
      if (gam <= x*(1.-x)) break;
      x = xmax;
    }
    Ene=x;
    if ( Ene >= (1.-Ee)) break;
    Ene = 1.-Ee;
  }
  G4double Enm=(2.-Ee-Ene);

  // initialisation of rotation parameters

  G4double costheta,sintheta,rphi,rtheta,rpsi;
  costheta= 1.-2./Ee-2./Ene+2./Ene/Ee;
  sintheta=std::sqrt(1.-costheta*costheta);
  

  rphi=twopi*G4UniformRand()*rad;
  rtheta=(std::acos(2.*G4UniformRand()-1.));
  rpsi=twopi*G4UniformRand()*rad;

  G4RotationMatrix rot;
  rot.set(rphi,rtheta,rpsi);

  // electron 0
  daughtermomentum[0]=std::sqrt(Ee*Ee*EMax*EMax+2.0*Ee*EMax * daughtermass[0]);
  G4ThreeVector direction0(0.0,0.0,1.0);

  direction0 *= rot;

  G4DynamicParticle* daughterparticle
    = new G4DynamicParticle (G4MT_daughters[0], direction0*daughtermomentum[0]);

  products->PushProducts(daughterparticle);
  
  // electronic neutrino  1

  daughtermomentum[1] = std::sqrt(Ene*Ene*EMax*EMax
                                  +2.0*Ene*EMax*daughtermass[1]);
  G4ThreeVector direction1(sintheta, 0.0, costheta);

  direction1 *= rot;

  G4DynamicParticle* daughterparticle1
    = new G4DynamicParticle (G4MT_daughters[1], direction1*daughtermomentum[1]);
  products->PushProducts(daughterparticle1);

  // muonnic neutrino 2
  
  daughtermomentum[2] = std::sqrt(Enm*Enm*EMax*EMax
                                  +2.0*Enm*EMax*daughtermass[2]);
  G4ThreeVector direction2(-Ene/Enm*sintheta,0,-Ee/Enm-Ene/Enm*costheta);

  direction2 *= rot;

  G4DynamicParticle* daughterparticle2
    = new G4DynamicParticle (G4MT_daughters[2], direction2*daughtermomentum[2]);
  products->PushProducts(daughterparticle2);

  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
  {
    G4cout << "G4MuonDecayChannel::DecayIt()";
    G4cout << " create decay products in rest frame " << G4endl;
    products->DumpInfo();
  }
#endif
  return products;
}
