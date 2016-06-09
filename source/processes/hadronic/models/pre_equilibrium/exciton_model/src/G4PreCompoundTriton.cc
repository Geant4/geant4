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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundTriton
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada add choice of options  
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
//
 
#include "G4PreCompoundTriton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Triton.hh"

G4PreCompoundTriton::G4PreCompoundTriton()
  : G4PreCompoundIon(G4Triton::Triton(), &theTritonCoulombBarrier)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentAthrd = ResidualAthrd;
  FragmentA = theA + ResidualA;
}

G4PreCompoundTriton::~G4PreCompoundTriton()
{}

G4double G4PreCompoundTriton::FactorialFactor(G4int N, const G4int P)
{
  return G4double((N-3)*(P-2)*(N-2)*(P-1)*(N-1)*P)/6.0; 
}
  
G4double G4PreCompoundTriton::CoalescenceFactor(G4int A)
{
  return 243.0/G4double(A*A);
}    

G4double G4PreCompoundTriton::GetRj(G4int nParticles, G4int nCharged)
{
  G4double rj = 0.0;
  if(nCharged >= 1 && (nParticles-nCharged) >= 2) {
    G4double denominator = 
      G4double(nParticles*(nParticles-1)*(nParticles-2));
    rj = G4double(3*nCharged*(nParticles-nCharged)*(nParticles-nCharged-1))
      /denominator; 
  }
  return rj;
}

//////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1,2 Chatterjee's paramaterization 
//OPT=3,4 Kalbach's parameterization 
// 
G4double G4PreCompoundTriton::CrossSection(G4double K)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentA = theA + ResidualA;
  FragmentAthrd = g4pow->Z13(FragmentA);

  if (OPTxs==0) { return GetOpt0( K); }
  else if( OPTxs==1 || OPTxs==2) { return GetOpt12( K); }
  else if (OPTxs==3 || OPTxs==4) { return GetOpt34( K); }
  else{
    std::ostringstream errOs;
    errOs << "BAD TRITON CROSS SECTION OPTION !!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;
  }
}

G4double G4PreCompoundTriton::GetAlpha()
{
  G4double C = 0.0;
  G4int aZ = theZ + ResidualZ;
  if (aZ >= 70) 
    {
      C = 0.10;
    } 
  else 
    {
      C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375; 
    }
 
  return 1.0 + C/3.0;
}

//
//********************* OPT=1,2 : Chatterjee's cross section *****************
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4PreCompoundTriton::GetOpt12(G4double K)
{
  G4double Kc=K;

  // JMQ xsec is set constat above limit of validity
  if (K > 50*MeV) { Kc=50*MeV; }

  G4double landa ,mu ,nu ,p , Ec,q,r,ji,xs;
 
  G4double    p0 = -11.04;
  G4double    p1 = 619.1;
  G4double    p2 = -2147.;
  G4double    landa0 = -0.0426;
  G4double    landa1 = -10.33;
  G4double    mm0 = 601.9;
  G4double    mu1 = 0.37;
  G4double    nu0 = 583.0;
  G4double    nu1 = -546.2;
  G4double    nu2 = 1.718;  
  G4double    delta=1.2;            

  Ec = 1.44*theZ*ResidualZ/(1.5*ResidualAthrd+delta);
  p = p0 + p1/Ec + p2/(Ec*Ec);
  landa = landa0*ResidualA + landa1;

  G4double resmu1 = g4pow->powZ(ResidualA,mu1); 
  mu = mm0*resmu1;
  nu = resmu1*(nu0 + nu1*Ec + nu2*(Ec*Ec));
  q = landa - nu/(Ec*Ec) - 2*p*Ec;
  r = mu + 2*nu/Ec + p*(Ec*Ec);
  
  ji=std::max(Kc,Ec);
  if(Kc < Ec) { xs = p*Kc*Kc + q*Kc + r;}
  else {xs = p*(Kc - ji)*(Kc - ji) + landa*Kc + mu + nu*(2 - Kc/ji)/ji ;}
                 
  if (xs <0.0) {xs=0.0;}
              
  return xs;
}

// *********** OPT=3,4 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundTriton::GetOpt34(G4double K)
//     ** t from o.m. of hafele, flynn et al
{
  G4double landa, mu, nu, p , signor(1.),sig;
  G4double ec,ecsq,xnulam,etest(0.),a; 
  G4double b,ecut,cut,ecut2,geom,elab;

  G4double     flow = 1.e-18;
  G4double     spill= 1.e+18;

  G4double     p0 = -21.45;
  G4double     p1 = 484.7;
  G4double     p2 = -1608.;
  G4double     landa0 = 0.0186;
  G4double     landa1 = -8.90;
  G4double     mm0 = 686.3;
  G4double     mu1 = 0.325;
  G4double     nu0 = 368.9;
  G4double     nu1 = -522.2;
  G4double     nu2 = -4.998;  
  
  G4double      ra=0.80;
        
  //JMQ 13/02/09 increase of reduced radius to lower the barrier
  // ec = 1.44 * theZ * ResidualZ / (1.5*ResidualAthrd+ra);
  ec = 1.44 * theZ * ResidualZ / (1.7*ResidualAthrd+ra);
  ecsq = ec * ec;
  p = p0 + p1/ec + p2/ecsq;
  landa = landa0*ResidualA + landa1;
  a = g4pow->powZ(ResidualA,mu1);
  mu = mm0 * a;
  nu = a* (nu0+nu1*ec+nu2*ecsq);  
  xnulam = nu / landa;
  if (xnulam > spill) { xnulam=0.; }
  if (xnulam >= flow) { etest = 1.2 *std::sqrt(xnulam); }
 
  a = -2.*p*ec + landa - nu/ecsq;
  b = p*ecsq + mu + 2.*nu/ec;
  ecut = 0.;
  cut = a*a - 4.*p*b;
  if (cut > 0.) { ecut = std::sqrt(cut); }
  ecut = (ecut-a) / (p+p);
  ecut2 = ecut;
  //JMQ 290310 for avoiding unphysical increase below minimum (at ecut)
  // ecut<0 means that there is no cut with energy axis, i.e. xs is set 
  // to 0 bellow minimum
  //  if (cut < 0.) ecut2 = ecut - 2.;
  if (cut < 0.) { ecut2 = ecut; }
  elab = K * FragmentA / G4double(ResidualA);
  sig = 0.;
 
  if (elab <= ec) { //start for E<Ec
    if (elab > ecut2) { sig = (p*elab*elab+a*elab+b) * signor; }
  }           //end for E<Ec
  else {           //start for E>Ec
    sig = (landa*elab+mu+nu/elab) * signor;
    geom = 0.;
    if (xnulam < flow || elab < etest) { return sig; }
    geom = std::sqrt(theA*K);
    geom = 1.23*ResidualAthrd + ra + 4.573/geom;
    geom = 31.416 * geom * geom;
    sig = std::max(geom,sig);
  }           //end for E>Ec
  return sig;
}
