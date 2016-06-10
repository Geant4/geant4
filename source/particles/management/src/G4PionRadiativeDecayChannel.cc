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
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History:
//               01 August 2007 P.Gumplinger
//               Reference: TRIUMF PIENU Technote:
//                          M. Blecher - "Inclusion of pi->enug in MC "
//                              Rate is for gammas > 100keV
//
// ------------------------------------------------------------
//
//
//

#include "G4PionRadiativeDecayChannel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"

G4PionRadiativeDecayChannel::G4PionRadiativeDecayChannel()
  : G4VDecayChannel(),
    beta(0.),cib(0.),csdp(0.),csdm(0.),cif(0.),cig(0.),
    xl(0.), yl(0.), xu(0.), yu(0.), d2wmax(0.)
{
}

G4PionRadiativeDecayChannel::
           G4PionRadiativeDecayChannel(const G4String& theParentName,
                                       G4double        theBR)
                            : G4VDecayChannel("Radiative Pion Decay",1)
{
  // set names for daughter particles
  if (theParentName == "pi+") {
    SetBR(theBR);
    SetParent("pi+");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "gamma");
    SetDaughter(2, "nu_e");
  } else if (theParentName == "pi-") {
    SetBR(theBR);
    SetParent("pi-");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "gamma");
    SetDaughter(2, "anti_nu_e");
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4RadiativePionDecayChannel:: constructor :";
      G4cout << " parent particle is not charged pion but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }

  beta = 3.6612e-03;

  cib  = 1.16141e-03;
  csdp = 3.45055e-02;
  csdm = 5.14122e-03;
  cif  = 4.63543e-05;
  cig  = 1.78928e-05;

  xl = 2.*0.1*MeV/139.57*MeV;
  yl = ((1.-xl) + std::sqrt((1-xl)*(1-xl)+4*beta*beta))/2.;

  xu = 1. - (yl - std::sqrt(yl*yl-4.*beta*beta))/2.;
  yu = 1. + beta*beta;

  d2wmax = D2W(xl,yl);

}

G4PionRadiativeDecayChannel::~G4PionRadiativeDecayChannel()
{
}
G4PionRadiativeDecayChannel::G4PionRadiativeDecayChannel(const G4PionRadiativeDecayChannel &right)
  :G4VDecayChannel(right),
   beta(right.beta),cib(right.cib),csdp(right.csdp),
   csdm(right.csdm),cif(right.cif),cig(right.cig),
   xl(right.xl), yl(right.yl), xu(right.xu), yu(right.yu), 
   d2wmax(right.d2wmax)
{
}

G4PionRadiativeDecayChannel & G4PionRadiativeDecayChannel::operator=(const G4PionRadiativeDecayChannel & right)
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
    beta = right.beta;
    cib  = right.cib;
    csdp = right.csdp;
    csdm = right.csdm;
    cif  = right.cif;
    cig  = right.cig;
    xl   = right.xl;
    yl   = right.yl;
    xu   = right.xu;
    yu   = right.yu; 
    d2wmax = right.d2wmax;
  }
  return *this;
}

G4DecayProducts *G4PionRadiativeDecayChannel::DecayIt(G4double) 
{

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) 
                 G4cout << "G4PionRadiativeDecayChannel::DecayIt ";
#endif

  if (G4MT_parent == 0) FillParent();  
  if (G4MT_daughters == 0) FillDaughters();

  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();

  G4double EMPI = parentmass;

  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }

  G4double EMASS = daughtermass[0];

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = 
                               new G4DynamicParticle( G4MT_parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4double x, y, d2w;

  do {

     do {

        x = xl + G4UniformRand()*(xu-xl);
        y = yl + G4UniformRand()*(yu-yl);

     } while (x+y <= 1.);

     d2w = D2W(x,y);

  } while (d2w <= G4UniformRand()*d2wmax);

//-----------------------------------------------------------------------
//
//      Calculate the angle between positron and photon (cosine)
//
  G4double cthetaGE =  (y*(x-2.)+2.*(1.-x+beta*beta)) /
                       (x*std::sqrt(y*y-4.*beta*beta));

//
//-----------------------------------------------------------------------
//
  G4double G = x * EMPI/2.;
  G4double E = y * EMPI/2.;
//
//-----------------------------------------------------------------------
//

  if (E < EMASS) E = EMASS;

  // calculate daughter momentum
  G4double daughtermomentum[2];

  daughtermomentum[0] = std::sqrt(E*E - EMASS*EMASS);

  G4double cthetaE = 2.*G4UniformRand()-1.;
  G4double sthetaE = std::sqrt(1.-cthetaE*cthetaE);

  G4double phiE = twopi*G4UniformRand()*rad;
  G4double cphiE = std::cos(phiE);
  G4double sphiE = std::sin(phiE);

  //Coordinates of the decay positron

  G4double px = sthetaE*cphiE;
  G4double py = sthetaE*sphiE;
  G4double pz = cthetaE;

  G4ThreeVector direction0(px,py,pz);

  G4DynamicParticle * daughterparticle0 
    = new G4DynamicParticle( G4MT_daughters[0], daughtermomentum[0]*direction0);

  products->PushProducts(daughterparticle0);

  daughtermomentum[1] = G;

  G4double sthetaGE = std::sqrt(1.-cthetaGE*cthetaGE);

  G4double phiGE = twopi*G4UniformRand()*rad;
  G4double cphiGE = std::cos(phiGE);
  G4double sphiGE = std::sin(phiGE);

  //Coordinates of the decay gamma with respect to the decay positron

  px = sthetaGE*cphiGE;
  py = sthetaGE*sphiGE;
  pz = cthetaGE;

  G4ThreeVector direction1(px,py,pz);

  direction1.rotateUz(direction0);

  G4DynamicParticle * daughterparticle1
    = new G4DynamicParticle( G4MT_daughters[1], daughtermomentum[1]*direction1);

  products->PushProducts(daughterparticle1);

// output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4PionRadiativeDecayChannel::DecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    products->DumpInfo();
  }
#endif

  return products;

}
