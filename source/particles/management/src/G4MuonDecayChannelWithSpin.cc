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
// G4MuonDecayChannelWithSpin class implementation
//
// References:
// - Florian Scheck "Muon Physics", in Physics Reports
//   (Review Section of Physics Letters) 44, No. 4 (1978)
//   187-248. North-Holland Publishing Company, Amsterdam at page 210 cc.
// - W.E. Fisher and F. Scheck, Nucl. Phys. B83 (1974) 25.

// Authors: P.Gumplinger and T.MacPhail, 17 August 2004 
// --------------------------------------------------------------------

#include "G4MuonDecayChannelWithSpin.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"

G4MuonDecayChannelWithSpin::G4MuonDecayChannelWithSpin()
  : G4MuonDecayChannel()
{
}

G4MuonDecayChannelWithSpin::
G4MuonDecayChannelWithSpin(const G4String& theParentName, 
                                 G4double  theBR)
  : G4MuonDecayChannel(theParentName,theBR)
{
}

G4MuonDecayChannelWithSpin::~G4MuonDecayChannelWithSpin()
{
}

G4MuonDecayChannelWithSpin::
G4MuonDecayChannelWithSpin(const G4MuonDecayChannelWithSpin& right)
  : G4MuonDecayChannel(right)
{
}

G4MuonDecayChannelWithSpin& G4MuonDecayChannelWithSpin::
operator=(const G4MuonDecayChannelWithSpin& right)
{
  if (this != &right)
  { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    delete parent_name;
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters > 0 )
    {
      daughters_name = new G4String*[numberOfDaughters];
      // copy daughters name
      for (G4int index=0; index<numberOfDaughters; ++index)
      {
        daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }
  return *this;
}

G4DecayProducts* G4MuonDecayChannelWithSpin::DecayIt(G4double) 
{
  // This version assumes V-A coupling with 1st order radiative correctons,
  //              the standard model Michel parameter values, but
  //              gives incorrect energy spectrum for neutrinos

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4MuonDecayChannelWithSpin::DecayIt ";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();

  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();

  G4double EMMU = parentmass;

  //daughters'mass
  G4double daughtermass[3]; 
  //G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; ++index)
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

  // calculate electron energy

  G4double michel_rho   = 0.75; //Standard Model Michel rho
  G4double michel_delta = 0.75; //Standard Model Michel delta
  G4double michel_xsi   = 1.00; //Standard Model Michel xsi
  G4double michel_eta   = 0.00; //Standard Model eta

  G4double rndm, x, ctheta;

  G4double FG; 
  G4double FG_max = 2.00;

  G4double W_mue  = (EMMU*EMMU+EMASS*EMASS)/(2.*EMMU);
  G4double x0     =           EMASS/W_mue;

  G4double x0_squared = x0*x0;
  
  // ***************************************************
  //     x0 <= x <= 1.   and   -1 <= y <= 1
  //
  //     F(x,y) = f(x)*g(x,y);   g(x,y) = 1.+g(x)*y
  // ***************************************************

  // ***** sampling F(x,y) directly (brute force) *****

  const std::size_t MAX_LOOP=10000;
  for (std::size_t loop_count=0; loop_count<MAX_LOOP; ++loop_count)
  {
    // Sample the positron energy by sampling from F

    rndm = G4UniformRand();

    x = x0 + rndm*(1.-x0);

    G4double x_squared = x*x;

    G4double F_IS, F_AS, G_IS, G_AS;

    F_IS = 1./6.*(-2.*x_squared+3.*x-x0_squared);
    F_AS = 1./6.*std::sqrt(x_squared-x0_squared)
                *(2.*x-2.+std::sqrt(1.-x0_squared));

    G_IS = 2./9.*(michel_rho-0.75)*(4.*x_squared-3.*x-x0_squared);
    G_IS = G_IS + michel_eta*(1.-x)*x0;

    G_AS = 3.*(michel_xsi-1.)*(1.-x);
    G_AS = G_AS+2.*(michel_xsi*michel_delta-0.75)
                  *(4.*x-4.+std::sqrt(1.-x0_squared));
    G_AS = 1./9.*std::sqrt(x_squared-x0_squared)*G_AS;

    F_IS = F_IS + G_IS;
    F_AS = F_AS + G_AS;

    // *** Radiative Corrections ***

    const G4double omega =  std::log(EMMU/EMASS);
    G4double R_IS = F_c(x,x0,omega);

    G4double F = 6.*F_IS + R_IS/std::sqrt(x_squared-x0_squared);

    // *** Radiative Corrections ***

    G4double R_AS = F_theta(x,x0,omega);

    rndm = G4UniformRand();

    ctheta = 2.*rndm-1.;

    G4double G = 6.*F_AS - R_AS/std::sqrt(x_squared-x0_squared);

    FG = std::sqrt(x_squared-x0_squared)*F*(1.+(G/F)*ctheta);

    if(FG>FG_max)
    {
      G4Exception("G4MuonDecayChannelWithSpin::DecayIt()",
                  "PART113", JustWarning,
                  "Problem in Muon Decay: FG > FG_max");
      FG_max = FG;
    }

    rndm = G4UniformRand();

    if (FG >= rndm*FG_max) break;
  }

  G4double energy = x * W_mue;

  rndm = G4UniformRand();

  G4double phi = twopi * rndm;

  if(energy < EMASS) energy = EMASS;

  // Calculate daughter momentum
  G4double daughtermomentum[3];

  daughtermomentum[0] = std::sqrt(energy*energy - EMASS*EMASS);

  G4double stheta = std::sqrt(1.-ctheta*ctheta);
  G4double cphi = std::cos(phi);
  G4double sphi = std::sin(phi);

  // Coordinates of the decay positron with respect to the muon spin
  G4double px = stheta*cphi;
  G4double py = stheta*sphi;
  G4double pz = ctheta;

  G4ThreeVector direction0(px,py,pz);

  direction0.rotateUz(parent_polarization);

  G4DynamicParticle * daughterparticle0 
    = new G4DynamicParticle( G4MT_daughters[0], daughtermomentum[0]*direction0);

  products->PushProducts(daughterparticle0);


  // daughter 1 ,2 (neutrinos)
  // create neutrinos in the C.M frame of two neutrinos
  G4double energy2 = parentmass-energy; 
  G4double vmass = std::sqrt((energy2-daughtermomentum[0])
                 * (energy2+daughtermomentum[0]));
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
  if (GetVerboseLevel()>1)
  {
    G4cout << "G4MuonDecayChannelWithSpin::DecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    G4double TT = daughterparticle0->GetTotalEnergy()
      + daughterparticle1->GetTotalEnergy()
      + daughterparticle2->GetTotalEnergy(); 
    G4cout << "e  " << daughterparticle0->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "nu1" << daughterparticle1->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "nu2" << daughterparticle2->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "total" << (TT-parentmass)/keV << G4endl;
    if (GetVerboseLevel()>2) { products->DumpInfo(); }
  }
#endif

  return products;
}

G4double G4MuonDecayChannelWithSpin::R_c(G4double x,G4double omega)
{
  G4int n_max = (G4int)(100.*x);

  if(n_max<10)n_max=10;

  G4double L2 = 0.0;

  for(G4int n=1; n<=n_max; ++n)
  {
    L2 += std::pow(x,n)/(n*n);
  }

  G4double r_c;

  r_c = 2.*L2-(pi*pi/3.)-2.;
  r_c = r_c + omega * (1.5+2.*std::log((1.-x)/x));
  r_c = r_c - std::log(x)*(2.*std::log(x)-1.);
  r_c = r_c + (3.*std::log(x)-1.-1./x)*std::log(1.-x);

  return r_c;
}
