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
// G4PionRadiativeDecayChannel class implementation
//      GEANT 4 class header file
//
// Author: P.Gumplinger, 30 July 2007 
// Reference: M. Blecher, TRIUMF/PIENU Technote
//            "Inclusion of pi->enug in the Monte Carlo" 
// --------------------------------------------------------------------

#include "G4PionRadiativeDecayChannel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"

namespace
{
  const G4double beta = 3.6612e-03;
  const G4double cib  = 1.16141e-03;
  const G4double csdp = 3.45055e-02;
  const G4double csdm = 5.14122e-03;
  const G4double cif  = 4.63543e-05;
  const G4double cig  = 1.78928e-05;
  const G4double xl = 2.*0.1*MeV/139.57*MeV;
  const G4double yl = ((1.-xl) + std::sqrt((1-xl)*(1-xl)+4*beta*beta))/2.;

  const G4double xu = 1. - (yl - std::sqrt(yl*yl-4.*beta*beta))/2.;
  const G4double yu = 1. + beta*beta;

  inline G4double D2W(const G4double x,const G4double y)
  {
    return cib*(1.-y)*(1.+((1.-x)*(1.-x)))/((x*x)*(x+y-1.)) +
        csdp*(1.-x)*((x+y-1.)*(x+y-1.)) +
        csdm*(1.-x)*((1.-y)*(1.-y)) +
        cif*(x-1.)*(1.-y)/x +
        cig*(1.-y)*(1.-x+(x*x)/(x+y-1.))/x;
  }

  const G4double d2wmax = D2W(xl,yl);
}

// --------------------------------------------------------------------
G4PionRadiativeDecayChannel::G4PionRadiativeDecayChannel()
  : G4VDecayChannel()
{
}

// --------------------------------------------------------------------
G4PionRadiativeDecayChannel::
G4PionRadiativeDecayChannel(const G4String& theParentName,
                                  G4double        theBR)
  : G4VDecayChannel("Radiative Pion Decay", 1)
{
  // set names for daughter particles
  if (theParentName == "pi+")
  {
    SetBR(theBR);
    SetParent("pi+");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "gamma");
    SetDaughter(2, "nu_e");
  }
  else if (theParentName == "pi-")
  {
    SetBR(theBR);
    SetParent("pi-");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "gamma");
    SetDaughter(2, "anti_nu_e");
  }
  else
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4RadiativePionDecayChannel::G4PionRadiativeDecayChannel()"
             << G4endl;
      G4cout << "Parent particle is not charged pion: ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

// --------------------------------------------------------------------
G4PionRadiativeDecayChannel::~G4PionRadiativeDecayChannel()
{
}

// --------------------------------------------------------------------
G4PionRadiativeDecayChannel::
G4PionRadiativeDecayChannel(const G4PionRadiativeDecayChannel& right)
  : G4VDecayChannel(right)
{
}

G4PionRadiativeDecayChannel&
G4PionRadiativeDecayChannel::operator=(const G4PionRadiativeDecayChannel& right)
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
      if (daughters_name != nullptr) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      //copy daughters name
      for (G4int index=0; index<numberOfDaughters; ++index)
      {
        daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }
  return *this;
}

// --------------------------------------------------------------------
G4DecayProducts* G4PionRadiativeDecayChannel::DecayIt(G4double) 
{

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) 
     G4cout << "G4PionRadiativeDecayChannel::DecayIt ";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();

  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();

  G4double EMPI = parentmass;

  // daughters'mass
  const G4int N_DAUGHTER=3;
  G4double daughtermass[N_DAUGHTER]; 
  //G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<N_DAUGHTER; ++index)
  {
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    //sumofdaughtermass += daughtermass[index];
  }

  G4double EMASS = daughtermass[0];

  // create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle* parentparticle
    = new G4DynamicParticle( G4MT_parent, dummy, 0.0);
  // create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4double x, y;

  const std::size_t MAX_LOOP=1000;

  for (std::size_t loop_counter1=0; loop_counter1<MAX_LOOP; ++loop_counter1)
  {
    for (std::size_t loop_counter2=0; loop_counter2<MAX_LOOP; ++loop_counter2)
    {
      x = xl + G4UniformRand()*(xu-xl);
      y = yl + G4UniformRand()*(yu-yl);
      if (x+y > 1.) break;
    }
    G4double d2w = D2W(x,y);
    if (d2w > G4UniformRand()*d2wmax) break;
  }

  // Calculate the angle between positron and photon (cosine)
  //
  G4double cthetaGE = (y*(x-2.)+2.*(1.-x+beta*beta)) /
                      (x*std::sqrt(y*y-4.*beta*beta));

  G4double G = x * EMPI/2.;
  G4double E = y * EMPI/2.;

  if (E < EMASS) E = EMASS;

  // calculate daughter momentum
  G4double daughtermomentum[2];

  daughtermomentum[0] = std::sqrt(E*E - EMASS*EMASS);

  G4double cthetaE = 2.*G4UniformRand()-1.;
  G4double sthetaE = std::sqrt(1.-cthetaE*cthetaE);

  G4double phiE = twopi*G4UniformRand()*rad;
  G4double cphiE = std::cos(phiE);
  G4double sphiE = std::sin(phiE);

  // Coordinates of the decay positron
  //
  G4double px = sthetaE*cphiE;
  G4double py = sthetaE*sphiE;
  G4double pz = cthetaE;

  G4ThreeVector direction0(px,py,pz);

  G4DynamicParticle * daughterparticle0 
    = new G4DynamicParticle(G4MT_daughters[0], daughtermomentum[0]*direction0);

  products->PushProducts(daughterparticle0);

  daughtermomentum[1] = G;

  G4double sthetaGE = std::sqrt(1.-cthetaGE*cthetaGE);

  G4double phiGE = twopi*G4UniformRand()*rad;
  G4double cphiGE = std::cos(phiGE);
  G4double sphiGE = std::sin(phiGE);

  // Coordinates of the decay gamma with respect to the decay positron
  //
  px = sthetaGE*cphiGE;
  py = sthetaGE*sphiGE;
  pz = cthetaGE;

  G4ThreeVector direction1(px,py,pz);

  direction1.rotateUz(direction0);

  G4DynamicParticle * daughterparticle1
    = new G4DynamicParticle(G4MT_daughters[1], daughtermomentum[1]*direction1);

  products->PushProducts(daughterparticle1);

  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
  {
    G4cout << "G4PionRadiativeDecayChannel::DecayIt() -";
    G4cout << " create decay products in rest frame " << G4endl;
    products->DumpInfo();
  }
#endif

  return products;
}
