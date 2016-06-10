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
// $Id: G4PreCompoundAlpha.cc 90591 2015-06-04 13:45:29Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundAlpha
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada add choice of options  
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
// 05.07.2013 J.M. Quesada FactorialFactor fixed
//

#include "G4PreCompoundAlpha.hh"
#include "G4SystemOfUnits.hh"
#include "G4Alpha.hh"

const G4double  p0 = 10.95;
const G4double  p1 = -85.2;
const G4double  p2 = 1146.;
const G4double  landa0 = 0.0643;
const G4double  landa1 = -13.96;
const G4double  mm0 = 781.2;
const G4double  mu1 = 0.29;
const G4double  nu0 = -304.7;
const G4double  nu1 = -470.0;
const G4double  nu2 = -8.580;        

G4PreCompoundAlpha::G4PreCompoundAlpha()
  : G4PreCompoundIon(G4Alpha::Alpha(), &theAlphaCoulombBarrier)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentAthrd = ResidualAthrd;
  FragmentA = theA + ResidualA;
}

G4PreCompoundAlpha::~G4PreCompoundAlpha()
{}

G4double G4PreCompoundAlpha::FactorialFactor(G4int N, G4int P)
{
  return G4double((N-4)*(P-3)*(N-3)*(P-2))*G4double((N-2)*(P-1)*(N-1)*P)/144.0;
}
  
G4double G4PreCompoundAlpha::CoalescenceFactor(G4int A)
{
  return 4096.0/G4double(A*A*A);  
}    

G4double G4PreCompoundAlpha::GetRj(G4int nParticles, G4int nCharged)
{
  G4double rj = 0.0;
  if(nCharged >=2 && (nParticles-nCharged) >=2 ) {
    G4double denominator = 
      G4double(nParticles*(nParticles-1)*(nParticles-2)*(nParticles-3));
    rj = 6.0*nCharged*(nCharged-1)*(nParticles-nCharged)*(nParticles-nCharged-1)
      /denominator;  
  }
  return rj;
}

/////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1,2 Chatterjee's paramaterization 
//OPT=3,4 Kalbach's parameterization 
// 
G4double G4PreCompoundAlpha::CrossSection(G4double K)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentA = theA + ResidualA;
  FragmentAthrd = g4pow->Z13(FragmentA);

  if (OPTxs==0)       { return GetOpt0( K); }
  else if(OPTxs <= 2) { return GetOpt12( K); }
  else                { return GetOpt34( K); }
}

G4double G4PreCompoundAlpha::GetAlpha()
{
  G4double C = 0.0;
  G4int aZ = theZ + ResidualZ;
  if (aZ <= 30) 
    {
      C = 0.10;
    } 
  else if (aZ <= 50) 
    {
      C = 0.1 - (aZ-30)*0.001;
    } 
  else if (aZ < 70) 
    {
      C = 0.08 - (aZ-50)*0.001;
    }
  else 
    {
      C = 0.06;
    }
  return 1.0+C;
}

//
//********************* OPT=1,2 : Chatterjee's cross section ********************
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4PreCompoundAlpha::GetOpt12(G4double K)
{
  // JMQ xsec is set constant above limit of validity
  G4double Kc = std::min(K, 50*MeV);
  G4double landa ,mu ,nu ,p , Ec,q,r,ji,xs;

  const G4double     delta=1.2;          

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
  
  xs = std::max(xs, 0.0);
  return xs;
}

// *********** OPT=3,4 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundAlpha::GetOpt34(G4double K)
// c     ** alpha from huizenga and igo
{
  const G4double     flow = 1.e-18;
  const G4double     spill= 1.e+18;  
  const G4double     ra = 1.20;
  const G4double     signor = 1.0;
        
  //JMQ 13/02/09 increase of reduced radius to lower the barrier
  // ec = 1.44 * theZ * ResidualZ / (1.5*ResidualAthrd+ra);
  G4double ec = 1.44 * theZ * ResidualZ / (1.7*ResidualAthrd+ra);
  G4double ecsq = ec * ec;
  G4double p = p0 + p1/ec + p2/ecsq;
  G4double landa = landa0*ResidualA + landa1;
  G4double a = g4pow->powZ(ResidualA,mu1);
  G4double mu = mm0 * a;
  G4double nu = a* (nu0+nu1*ec+nu2*ecsq);  
  G4double xnulam = nu / landa;
  G4double etest = 0.0;
  if (xnulam > spill)      { xnulam=0.; }
  else if (xnulam >= flow) { etest = 1.2 *std::sqrt(xnulam); }

  a = -2.*p*ec + landa - nu/ecsq;
  G4double b = p*ecsq + mu + 2.*nu/ec;
  G4double ecut = 0.;
  G4double cut = a*a - 4.*p*b;
  if (cut > 0.) { ecut = std::sqrt(cut); }
  ecut = (ecut-a) / (2*p);

  //JMQ 290310 for avoiding unphysical increase below minimum (at ecut)
  // ecut<0 means that there is no cut with energy axis, i.e. xs is set 
  // to 0 bellow minimum

  G4double elab = K * FragmentA / G4double(ResidualA);
  G4double sig = 0.;
  
  if (elab <= ec) { 
    if (elab > ecut) { sig = std::max(0.0,(p*elab*elab+a*elab+b) * signor); }

  } else {  
    sig = (landa*elab+mu+nu/elab) * signor;
    G4double geom = 0.;
    if (xnulam >= flow && elab >= etest) {
      geom = std::sqrt(theA*K);
      geom = 1.23*ResidualAthrd + ra + 4.573/geom;
      geom = 31.416 * geom * geom;
    }
    sig = std::max(sig, geom);
  }          
  return sig;
}
