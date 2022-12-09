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
// G4MuonRadiativeDecayChannelWithSpin class implementation
//
// References:
// - TRIUMF/TWIST Technote TN-55:
//   "Radiative muon decay" by P. Depommier and A. Vacheret
// - Yoshitaka Kuno and Yasuhiro Okada
//   "Muon Decays and Physics Beyond the Standard Model"
//   Rev. Mod. Phys. 73, 151 (2001)
//
// Author: P.Gumplinger - Triumf, 25 July 2007   
// Revision: D.Mingming - Center for HEP, Tsinghua Univ., 10 August 2011 
// --------------------------------------------------------------------

#include "G4MuonRadiativeDecayChannelWithSpin.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"

G4MuonRadiativeDecayChannelWithSpin::G4MuonRadiativeDecayChannelWithSpin()
  : G4VDecayChannel()
{
}

G4MuonRadiativeDecayChannelWithSpin::
G4MuonRadiativeDecayChannelWithSpin(const G4String& theParentName,
                                          G4double  theBR)
  : G4VDecayChannel("Radiative Muon Decay",1)
{
  // set names for daughter particles
  if (theParentName == "mu+")
  {
    SetBR(theBR);
    SetParent("mu+");
    SetNumberOfDaughters(4);
    SetDaughter(0, "e+");
    SetDaughter(1, "gamma");
    SetDaughter(2, "nu_e");
    SetDaughter(3, "anti_nu_mu");
  }
  else if (theParentName == "mu-")
  {
    SetBR(theBR);
    SetParent("mu-");
    SetNumberOfDaughters(4);
    SetDaughter(0, "e-");
    SetDaughter(1, "gamma");
    SetDaughter(2, "anti_nu_e");
    SetDaughter(3, "nu_mu");
  }
  else
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4RadiativeMuonDecayChannel::G4RadiativeMuonDecayChannel():";
      G4cout << " parent particle is not muon but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

G4MuonRadiativeDecayChannelWithSpin::~G4MuonRadiativeDecayChannelWithSpin()
{
}

G4MuonRadiativeDecayChannelWithSpin::
G4MuonRadiativeDecayChannelWithSpin(const G4MuonRadiativeDecayChannelWithSpin& r)
  : G4VDecayChannel(r)
{
}

G4MuonRadiativeDecayChannelWithSpin& G4MuonRadiativeDecayChannelWithSpin::
operator=(const G4MuonRadiativeDecayChannelWithSpin& right)
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
    if ( numberOfDaughters > 0 )
    {
      if (daughters_name != nullptr) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      // copy daughters name
      for (G4int index=0; index<numberOfDaughters; ++index)
      {
        daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
    parent_polarization = right.parent_polarization;
  }
  return *this;
}


G4DecayProducts* G4MuonRadiativeDecayChannelWithSpin::DecayIt(G4double) 
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) 
    G4cout << "G4MuonRadiativeDecayChannelWithSpin::DecayIt()";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();

  // parent mass
  G4double parentmass = G4MT_parent->GetPDGMass();

  G4double EMMU = parentmass;

  // daughters'mass
  G4double daughtermass[4]; 
  //G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<4; ++index)
  {
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    //sumofdaughtermass += daughtermass[index];
  }

  G4double EMASS = daughtermass[0];

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle* parentparticle
    = new G4DynamicParticle( G4MT_parent, dummy, 0.0 );

  // create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4double eps = EMASS/EMMU;

  G4double som0, x, y, xx, yy, zz;
  G4double cthetaE, cthetaG, cthetaGE, phiE, phiG;
  const std::size_t MAX_LOOP=10000;

  for (std::size_t loop_counter1=0; loop_counter1<MAX_LOOP; ++loop_counter1)
  {
    for (std::size_t loop_counter2=0; loop_counter2<MAX_LOOP; ++loop_counter2)
    {
      // -------------------------------------------------------------------
      // Build two vectors of random length and random direction, for the
      // positron and the photon.
      // x/y is the length of the vector, xx, yy and zz the components,
      // phi is the azimutal angle, theta the polar angle.
      // -------------------------------------------------------------------

      // For the positron
      //
      x = G4UniformRand();

      rn3dim(xx,yy,zz,x);

      if(std::fabs((xx*xx)+(yy*yy)+(zz*zz)-(x*x))>0.001)
      {
        G4cout << "Norm of x not correct" << G4endl;
      }

      phiE = atan4(xx,yy);
      cthetaE = zz/x;
      G4double sthetaE = std::sqrt((xx*xx)+(yy*yy))/x;

      // What you get:
      //
      // x       = positron energy
      // phiE    = azimutal angle of positron momentum
      // cthetaE = cosine of polar angle of positron momentum
      // sthetaE = sine of polar angle of positron momentum
      //
      //// G4cout << " x, xx, yy, zz " << x  << " " << xx << " " 
      ////                             << yy << " " << zz << G4endl;
      //// G4cout << " phiE, cthetaE, sthetaE " << phiE    << " "
      ////                                      << cthetaE << " " 
      ////                                      << sthetaE << " " << G4endl;

      // For the photon
      //
      y = G4UniformRand();

      rn3dim(xx,yy,zz,y);

      if(std::fabs((xx*xx)+(yy*yy)+(zz*zz)-(y*y))>0.001)
      {
        G4cout << " Norm of y not correct " << G4endl;
      }

      phiG = atan4(xx,yy);
      cthetaG = zz/y;
      G4double sthetaG = std::sqrt((xx*xx)+(yy*yy))/y;

      // What you get:
      //
      // y       = photon energy
      // phiG    = azimutal angle of photon momentum
      // cthetaG = cosine of polar angle of photon momentum
      // sthetaG = sine of polar angle of photon momentum
      //
      //// G4cout << " y, xx, yy, zz " << y  << " " << xx << " "
      ////                             << yy << " " << zz << G4endl;
      //// G4cout << " phiG, cthetaG, sthetaG " << phiG    << " "
      ////                                      << cthetaG << " "
      ////                                      << sthetaG << " " << G4endl;

      //      Calculate the angle between positron and photon (cosine)
      //
      cthetaGE = cthetaE*cthetaG+sthetaE*sthetaG*std::cos(phiE-phiG);

      //// G4cout << x << " " << cthetaE << " " << sthetaE << " "
      ////        << y << " " << cthetaG << " " << sthetaG << " "
      ////        << cthetaGE

      G4double term0 = eps*eps;
      G4double term1 = x*((1.0-eps)*(1.0-eps))+2.0*eps;
      G4double beta  = std::sqrt( x*((1.0-eps)*(1.0-eps))*
                                 (x*((1.0-eps)*(1.0-eps))+4.0*eps))/term1;
      G4double delta = 1.0-beta*cthetaGE;

      G4double term3 = y*(1.0-(eps*eps));
      G4double term6 = term1*delta*term3;

      G4double Qsqr = (1.0-term1-term3+term0+0.5*term6)/((1.0-eps)*(1.0-eps));

      // Check the kinematics.
      //
      if  ( Qsqr>=0.0 && Qsqr<=1.0 ) break;

    } // end loop count

    // Do the calculation for -1 muon polarization (i.e. mu+)
    //
    G4double Pmu = -1.0;
    if (GetParentName() == "mu-")  { Pmu = +1.0; }

    som0 = fron(Pmu,x,y,cthetaE,cthetaG,cthetaGE);

    // Sample the decay rate
    //
    if (G4UniformRand()*177.0 <= som0) break;
  }

  G4double E = EMMU/2.*(x*((1.-eps)*(1.-eps))+2.*eps);
  G4double G = EMMU/2.*y*(1.-eps*eps);

  if(E < EMASS) E = EMASS;

  // calculate daughter momentum
  G4double daughtermomentum[4];

  daughtermomentum[0] = std::sqrt(E*E - EMASS*EMASS);

  G4double sthetaE = std::sqrt(1.-cthetaE*cthetaE);
  G4double cphiE = std::cos(phiE);
  G4double sphiE = std::sin(phiE);

  // Coordinates of the decay positron with respect to the muon spin

  G4double px = sthetaE*cphiE;
  G4double py = sthetaE*sphiE;
  G4double pz = cthetaE;

  G4ThreeVector direction0(px,py,pz);

  direction0.rotateUz(parent_polarization);

  G4DynamicParticle * daughterparticle0 
    = new G4DynamicParticle( G4MT_daughters[0], daughtermomentum[0]*direction0);

  products->PushProducts(daughterparticle0);

  daughtermomentum[1] = G;

  G4double sthetaG = std::sqrt(1.-cthetaG*cthetaG);
  G4double cphiG = std::cos(phiG);
  G4double sphiG = std::sin(phiG);

  // Coordinates of the decay gamma with respect to the muon spin

  px = sthetaG*cphiG;
  py = sthetaG*sphiG;
  pz = cthetaG;

  G4ThreeVector direction1(px,py,pz);

  direction1.rotateUz(parent_polarization);

  G4DynamicParticle * daughterparticle1
    = new G4DynamicParticle( G4MT_daughters[1], daughtermomentum[1]*direction1);

  products->PushProducts(daughterparticle1);

  // daughter 3 ,4 (neutrinos)
  // create neutrinos in the C.M frame of two neutrinos

  G4double energy2 = parentmass-E-G;

  G4ThreeVector P34 = -1.*(daughtermomentum[0]*direction0
                          +daughtermomentum[1]*direction1);  
  G4double vmass2 = energy2*energy2 - P34.mag2();
  G4double vmass = std::sqrt(vmass2);

  G4double costhetan = 2.*G4UniformRand()-1.0;
  G4double sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
  G4double phin  = twopi*G4UniformRand()*rad;
  G4double sinphin = std::sin(phin);
  G4double cosphin = std::cos(phin);

  G4ThreeVector direction2(sinthetan*cosphin,sinthetan*sinphin,costhetan);

  G4DynamicParticle * daughterparticle2
    = new G4DynamicParticle( G4MT_daughters[2], direction2*(vmass/2.));
  G4DynamicParticle * daughterparticle3
    = new G4DynamicParticle( G4MT_daughters[3], direction2*(-1.0*vmass/2.));

  // boost to the muon rest frame
  G4double beta = P34.mag()/energy2;
  G4ThreeVector direction34 = P34.unit();

  G4LorentzVector p4 = daughterparticle2->Get4Momentum();
  p4.boost(direction34.x()*beta,direction34.y()*beta,direction34.z()*beta);
  daughterparticle2->Set4Momentum(p4);

  p4 = daughterparticle3->Get4Momentum();
  p4.boost(direction34.x()*beta,direction34.y()*beta,direction34.z()*beta);
  daughterparticle3->Set4Momentum(p4);

  products->PushProducts(daughterparticle2);
  products->PushProducts(daughterparticle3);

  daughtermomentum[2] = daughterparticle2->GetTotalMomentum();
  daughtermomentum[3] = daughterparticle3->GetTotalMomentum();

  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
  {
    G4cout << "G4MuonRadiativeDecayChannelWithSpin::DecayIt() -";
    G4cout << " create decay products in rest frame " <<G4endl;
    G4double TT = daughterparticle0->GetTotalEnergy()
      + daughterparticle1->GetTotalEnergy()
      + daughterparticle2->GetTotalEnergy() 
      + daughterparticle3->GetTotalEnergy(); 
    G4cout << "e    :" << daughterparticle0->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "gamma:" << daughterparticle1->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "nu2  :" << daughterparticle2->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "nu2  :" << daughterparticle3->GetTotalEnergy()/MeV << G4endl; 
    G4cout << "total:" << (TT-parentmass)/keV << G4endl;
    if (GetVerboseLevel()>1) { products->DumpInfo(); }
  }
#endif

 return products;
}

G4double G4MuonRadiativeDecayChannelWithSpin::fron(G4double Pmu,
                                                   G4double x,
                                                   G4double y,
                                                   G4double cthetaE,
                                                   G4double cthetaG,
                                                   G4double cthetaGE)
{
  G4double mu  = 105.65;
  G4double me  =   0.511;
  G4double rho =   0.75;
  G4double del =   0.75;
  G4double eps =   0.0;
  G4double kap =   0.0;
  G4double ksi =   1.0;

  G4double delta = 1-cthetaGE;

  // Calculation of the functions f(x,y)

  G4double f_1s  = 12.0*((y*y)*(1.0-y)+x*y*(2.0-3.0*y)
                   +2.0*(x*x)*(1.0-2.0*y)-2.0*(x*x*x));
  G4double f0s   = 6.0*(-x*y*(2.0-3.0*(y*y))
                   -2.0*(x*x)*(1.0-y-3.0*(y*y))+2.0*(x*x*x)*(1.0+2.0*y));
  G4double f1s   = 3.0*((x*x)*y*(2.0-3.0*y-3.0*(y*y))
                   -(x*x*x)*y*(4.0+3.0*y));
  G4double f2s   = 1.5*((x*x*x)*(y*y)*(2.0+y));

  G4double f_1se = 12.0*(x*y*(1.0-y)+(x*x)*(2.0-3.0*y)
                   -2.0*(x*x*x));
  G4double f0se  = 6.0*(-(x*x)*(2.0-y-2.0*(y*y))
                   +(x*x*x)*(2.0+3.0*y));
  G4double f1se  = -3.0*(x*x*x)*y*(2.0+y);
  G4double f2se  = 0.0;

  G4double f_1sg = 12.0*((y*y)*(1.0-y)+x*y*(1.0-2.0*y)
                   -(x*x)*y);
  G4double f0sg  = 6.0*(-x*(y*y)*(2.0-3.0*y)-(x*x)*y*(1.0-4.0*y)
                   +(x*x*x)*y);
  G4double f1sg  = 3.0*((x*x)*(y*y)*(1.0-3.0*y)
                   -2.0*(x*x*x)*(y*y));
  G4double f2sg  = 1.5*(x*x*x)*(y*y*y);

  G4double f_1v  = 8.0*((y*y)*(3.0-2.0*y)+6.0*x*y*(1.0-y)
                   +2.0*(x*x)*(3.0-4.0*y)-4.0*(x*x*x));
  G4double f0v   = 8.0*(-x*y*(3.0-y-(y*y))-(x*x)*(3.0-y-4.0*(y*y))
                   +2.0*(x*x*x)*(1.0+2.0*y));
  G4double f1v   = 2.0*((x*x)*y*(6.0-5.0*y-2.0*(y*y))
                   -2.0*(x*x*x)*y*(4.0+3.0*y));
  G4double f2v   = 2.0*(x*x*x)*(y*y)*(2.0+y);

  G4double f_1ve = 8.0*(x*y*(1.0-2.0*y)
                   +2.0*(x*x)*(1.0-3.0*y)-4.0*(x*x*x));
  G4double f0ve  = 4.0*(-(x*x)*(2.0-3.0*y-4.0*(y*y))
                   +2.0*(x*x*x)*(2.0+3.0*y));
  G4double f1ve  = -4.0*(x*x*x)*y*(2.0+y);
  G4double f2ve  = 0.0;

  G4double f_1vg = 8.0*((y*y)*(1.0-2.0*y)+x*y*(1.0-4.0*y)
                   -2.0*(x*x)*y);
  G4double f0vg  = 4.0*(2.0*x*(y*y)*(1.0+y)-(x*x)*y*(1.0-4.0*y)
                   +2.0*(x*x*x)*y);
  G4double f1vg  = 2.0*((x*x)*(y*y)*(1.0-2.0*y)
                   -4.0*(x*x*x)*(y*y));
  G4double f2vg  = 2.0*(x*x*x)*(y*y*y);

  G4double f_1t  = 8.0*((y*y)*(3.0-y)+3.0*x*y*(2.0-y)
                   +2.0*(x*x)*(3.0-2.0*y)-2.0*(x*x*x));
  G4double f0t   = 4.0*(-x*y*(6.0+(y*y))
                   -2.0*(x*x)*(3.0+y-3.0*(y*y))+2.0*(x*x*x)*(1.0+2.0*y));
  G4double f1t   = 2.0*((x*x)*y*(6.0-5.0*y+(y*y))
                   -(x*x*x)*y*(4.0+3.0*y));
  G4double f2t   = (x*x*x)*(y*y)*(2.0+y);

  G4double f_1te = -8.0*(x*y*(1.0+3.0*y)+(x*x)*(2.0+3.0*y)
                   +2.0*(x*x*x));
  G4double f0te  = 4.0*((x*x)*(2.0+3.0*y+4.0*(y*y))
                   +(x*x*x)*(2.0+3.0*y));
  G4double f1te  = -2.0*(x*x*x)*y*(2.0+y);
  G4double f2te  = 0.0;

  G4double f_1tg = -8.0*((y*y)*(1.0+y)+x*y+(x*x)*y);
  G4double f0tg  = 4.0*(x*(y*y)*(2.0-y)+(x*x)*y*(1.0+2.0*y)
                   +(x*x*x)*y);
  G4double f1tg  = -2.0*((x*x)*(y*y)*(1.0-y)+2.0*(x*x*x)*y);
  G4double f2tg  = (x*x*x)*(y*y*y);

  G4double term = delta+2.0*(me*me)/((mu*mu)*(x*x));
  term = 1.0/term;

  G4double nss = term*f_1s+f0s+delta*f1s+(delta*delta)*f2s;
  G4double nv = term*f_1v+f0v+delta*f1v+(delta*delta)*f2v;
  G4double nt = term*f_1t+f0t+delta*f1t+(delta*delta)*f2t;

  G4double nse = term*f_1se+f0se+delta*f1se+(delta*delta)*f2se;
  G4double nve = term*f_1ve+f0ve+delta*f1ve+(delta*delta)*f2ve;
  G4double nte = term*f_1te+f0te+delta*f1te+(delta*delta)*f2te;

  G4double nsg = term*f_1sg+f0sg+delta*f1sg+(delta*delta)*f2sg;
  G4double nvg = term*f_1vg+f0vg+delta*f1vg+(delta*delta)*f2vg;
  G4double ntg = term*f_1tg+f0tg+delta*f1tg+(delta*delta)*f2tg;

  G4double term1 = nv;
  G4double term2 = 2.0*nss+nv-nt;
  G4double term3 = 2.0*nss-2.0*nv+nt;

  G4double term1e = 1.0/3.0*(1.0-4.0/3.0*del);
  G4double term2e = 2.0*nse+5.0*nve-nte;
  G4double term3e = 2.0*nse-2.0*nve+nte;

  G4double term1g = 1.0/3.0*(1.0-4.0/3.0*del);
  G4double term2g = 2.0*nsg+5.0*nvg-ntg;
  G4double term3g = 2.0*nsg-2.0*nvg+ntg;

  G4double som00 = term1+(1.0-4.0/3.0*rho)*term2+eps*term3;
  G4double som01 = Pmu*ksi*(cthetaE*(nve-term1e*term2e+kap*term3e)
                   +cthetaG*(nvg-term1g*term2g+kap*term3g));

  G4double som0 = (som00+som01)/y;
  som0  = fine_structure_const/8./(twopi*twopi*twopi)*som0;

  return som0;
}
