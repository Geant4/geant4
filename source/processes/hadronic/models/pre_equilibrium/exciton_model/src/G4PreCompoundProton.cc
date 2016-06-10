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
// $Id: G4PreCompoundProton.cc 90591 2015-06-04 13:45:29Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundProton
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada added external choice of inverse cross section option
// 21.08.2008 J. M. Quesada added external choice for superimposed Coulomb 
//                          barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
//

#include "G4PreCompoundProton.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

G4PreCompoundProton::G4PreCompoundProton()
  : G4PreCompoundNucleon(G4Proton::Proton(), &theProtonCoulombBarrier)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentAthrd = ResidualAthrd;
  FragmentA = theA + ResidualA;
}

G4PreCompoundProton::~G4PreCompoundProton()
{}

G4double G4PreCompoundProton::GetRj(G4int nParticles, G4int nCharged)
{
  G4double rj = 0.0;
  if(nParticles > 0) { 
    rj = static_cast<G4double>(nCharged)/static_cast<G4double>(nParticles);
  }
  return rj;
}

////////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1 Chatterjee's paramaterization 
//OPT=2,4 Wellisch's parametarization
//OPT=3 Kalbach's parameterization 
// 
G4double G4PreCompoundProton::CrossSection(G4double K)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentA = theA + ResidualA;
  FragmentAthrd = g4pow->Z13(FragmentA);

  if (OPTxs==0)        { return GetOpt0( K); }
  else if( OPTxs == 1) { return GetOpt1( K); }
  else if( OPTxs == 2) { return GetOpt2( K); }
  else                 { return GetOpt3( K); }
}

G4double G4PreCompoundProton::GetAlpha()
{
  G4int aZ = ResidualZ;
  G4double C = 0.0;
  if (aZ >= 70) 
    {
      C = 0.10;
    } 
  else 
    {
      C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ 
	   - 0.66612e-01)*aZ + 0.98375;
    }
  return 1.0 + C;
}

G4double G4PreCompoundProton::GetBeta() 
{
  return -GetCoulombBarrier();
}
  
//********************* OPT=1 : Chatterjee's cross section *********************
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4PreCompoundProton::GetOpt1(G4double K)
{
  G4double Kc=K; 

  // JMQ  xsec is set constat above limit of validity
  if (K > 50*MeV) { Kc = 50*MeV; }

  const G4double p0 = 15.72;
  const G4double p1 = 9.65;
  const G4double p2 = -449.0;
  const G4double landa0 = 0.00437;
  const G4double landa1 = -16.58;
  const G4double mm0 = 244.7;
  const G4double mu1 = 0.503;
  const G4double nu0 = 273.1;
  const G4double nu1 = -182.4;
  const G4double nu2 = -1.872;  
  const G4double delta = 0.;  

  G4double Ec = 1.44*theZ*ResidualZ/(1.5*ResidualAthrd+delta);
  G4double p = p0 + p1/Ec + p2/(Ec*Ec);
  G4double landa = landa0*ResidualA + landa1;

  G4double resmu1 = g4pow->powZ(ResidualA,mu1); 
  G4double mu = mm0*resmu1;
  G4double nu = resmu1*(nu0 + nu1*Ec + nu2*(Ec*Ec));
  G4double q = landa - nu/(Ec*Ec) - 2*p*Ec;
  G4double r = mu + 2*nu/Ec + p*(Ec*Ec);

  G4double ji = std::max(Kc,Ec);
  G4double xs = 0.0;

  if(Kc < Ec) { xs = p*Kc*Kc + q*Kc + r;}
  else {xs = p*(Kc - ji)*(Kc - ji) + landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}

  xs = std::max(xs, 0.0);
  return xs; 
}

//************* OPT=2 : Welisch's proton reaction cross section ***************

G4double G4PreCompoundProton::GetOpt2(G4double K)
{
  // This is redundant when the Coulomb  barrier is overimposed to all 
  // cross sections 
  // It should be kept when Coulomb barrier only imposed at OPTxs=2

  if(!useSICB && K<=theCoulombBarrier) { return 0.0; }

  G4double eekin=K;
  G4int rnneu=ResidualA-ResidualZ;
  G4double ekin=eekin/1000;
  G4double r0=1.36*1.e-15;
  G4double fac=pi*r0*r0;
  G4double b0=2.247-0.915*(1.-1./ResidualAthrd);
  G4double fac1=b0*(1.-1./ResidualAthrd);
  G4double fac2=1.;
  if(rnneu > 1.5) { fac2 = g4pow->logZ(rnneu); }
  G4double xine_th= 1.e+31*fac*fac2*(1.+ResidualAthrd-fac1);
  xine_th=(1.-0.15*G4Exp(-ekin))*xine_th/(1.00-0.0007*ResidualA);	
  G4double ff1=0.70-0.0020*ResidualA;
  G4double ff2=1.00+1/G4double(ResidualA);
  G4double ff3=0.8+18/G4double(ResidualA)-0.002*ResidualA;
  G4double log10E = G4Log(ekin)/g4pow->logZ(10);
  fac=1.-(1./(1.+G4Exp(-8.*ff1*(log10E + 1.37*ff2))));
  xine_th=xine_th*(1.+ff3*fac);
  ff1=1.-1/G4double(ResidualA)-0.001*ResidualA;
  ff2=1.17-2.7/G4double(ResidualA)-0.0014*ResidualA;
  fac=-8.*ff1*(log10E + 2.0*ff2);
  xine_th /= (1.+G4Exp(fac));    
        
  xine_th = std::max(xine_th, 0.0);
  return xine_th;
}

// *********** OPT=3 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundProton::GetOpt3(const  G4double K)
{
  //     ** p from  becchetti and greenlees (but modified with sub-barrier
  //     ** correction function and xp2 changed from -449)

  const G4double p0 = 15.72;
  const G4double p1 = 9.65;
  const G4double p2 = -300.;
  const G4double landa0 = 0.00437;
  const G4double landa1 = -16.58;
  const G4double mm0 = 244.7;
  const G4double mu1 = 0.503;
  const G4double nu0 = 273.1;
  const G4double nu1 = -182.4;
  const G4double nu2 = -1.872;
  
  const G4double flow  = 1.e-18;
  const G4double spill = 1.e+18; 
  const G4double ra = 0.0; 
   
  G4double signor = 1.0;
  if (ResidualA <= 60)      { signor = 0.92; }
  else if (ResidualA < 100) { signor = 0.8 + ResidualA*0.002; }
  
  G4double ec = 1.44 * theZ * ResidualZ / (1.5*ResidualAthrd+ra);
  G4double ecsq = ec * ec;
  G4double p = p0 + p1/ec + p2/ecsq;
  G4double landa = landa0*ResidualA + landa1;
  G4double a = g4pow->powZ(ResidualA,mu1);
  G4double mu = mm0 * a;
  G4double nu = a* (nu0+nu1*ec+nu2*ecsq);
  
  G4double c =std::min(3.15,ec*0.5);
  G4double w = 0.7 * c / 3.15; 
  
  G4double etest = 0.0;
  G4double xnulam = nu / landa;
  if(xnulam > spill)      { xnulam=0.; }
  else if(xnulam >= flow) { etest = std::sqrt(xnulam) + 7.; }
  
  a = -2.*p*ec + landa - nu/ecsq;
  G4double b = p*ecsq + mu + 2.*nu/ec;
  G4double ecut = 0.;
  G4double cut = a*a - 4.*p*b;
  if (cut > 0.) { ecut = std::sqrt(cut); }
  ecut = (ecut-a) / (2*p);
 
  //JMQ 290310 for avoiding unphysical increase below minimum (at ecut)
  // ecut<0 means that there is no cut with energy axis, i.e. xs is set 
  // to 0 bellow minimum

  G4double elab = K * FragmentA /G4double(ResidualA);
  G4double sig = 0.;
  if (elab <= ec) { 
    if (elab > ecut) { sig = (p*elab*elab+a*elab+b) * signor; }
    
    G4double signor2 = (ec-elab-c) / w;
    sig /= (1. + G4Exp(signor2));

  } else { 
    sig = (landa*elab+mu+nu/elab) * signor;
    G4double geom = 0.;    
    if (xnulam >= flow && elab >= etest) { 
      geom = std::sqrt(theA*K);
      geom = 1.23*ResidualAthrd + ra + 4.573/geom;
      geom = 31.416 * geom * geom;
      sig = std::max(geom, sig);
    }
  }  
  sig = std::max(sig, 0.0);
  return sig;
}
