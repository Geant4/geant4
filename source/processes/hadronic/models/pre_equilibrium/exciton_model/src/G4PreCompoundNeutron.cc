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
// $Id: G4PreCompoundNeutron.cc 90591 2015-06-04 13:45:29Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundNeutron
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada add choice of options  
// 10.02.2009 J. M. Quesada set default opt3
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
// 

#include "G4PreCompoundNeutron.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"

G4PreCompoundNeutron::G4PreCompoundNeutron()
  : G4PreCompoundNucleon(G4Neutron::Neutron(), &theNeutronCoulombBarrier)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentAthrd = ResidualAthrd;
  FragmentA = theA + ResidualA;
}

G4PreCompoundNeutron::~G4PreCompoundNeutron()
{}

G4double G4PreCompoundNeutron::GetRj(G4int nParticles, G4int nCharged)
{
  G4double rj = 0.0;
  if(nParticles > 0) { 
    rj = static_cast<G4double>(nParticles - nCharged)/
      static_cast<G4double>(nParticles);
  }
  return rj;
}

////////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1,2 Chatterjee's paramaterization 
//OPT=3,4 Kalbach's parameterization 
// 
G4double G4PreCompoundNeutron::CrossSection(const  G4double K)
{
  ResidualA = GetRestA();
  ResidualZ = GetRestZ(); 
  theA = GetA();
  theZ = GetZ();
  ResidualAthrd = ResidualA13();
  FragmentA = theA + ResidualA;
  FragmentAthrd = g4pow->Z13(FragmentA);

  if (OPTxs==0)        { return GetOpt0( K); }
  else if( OPTxs <= 2) { return GetOpt12( K); }
  else                 { return GetOpt34( K); }
}

G4double G4PreCompoundNeutron::GetAlpha()
{
  return 0.76+2.2/ResidualAthrd;
}

G4double G4PreCompoundNeutron::GetBeta() 
{
  //   return (2.12/std::pow(GetRestA(),2.0/3.0)-0.05)*MeV/GetAlpha();
  return (2.12/(ResidualAthrd*ResidualAthrd)-0.05)*MeV/GetAlpha();
}

//********************* OPT=1,2 : Chatterjee's cross section ***************
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4PreCompoundNeutron::GetOpt12(G4double K)
{
  G4double Kc=K;

  // Pramana (Bechetti & Greenles) for neutrons is chosen 

  // JMQ  xsec is set constat above limit of validity
  if (K > 50*MeV) { Kc = 50*MeV; }

  const G4double landa0 = 18.57;
  const G4double landa1 = -22.93;
  const G4double mm0 = 381.7;
  const G4double mu1 = 24.31;
  const G4double nu0 = 0.172;
  const G4double nu1 = -15.39;
  const G4double nu2 = 804.8;

  G4double landa = landa0/ResidualAthrd + landa1;
  G4double mu = mm0*ResidualAthrd + mu1*ResidualAthrd*ResidualAthrd;
  G4double nu = nu0*ResidualAthrd*ResidualA + nu1*ResidualAthrd*ResidualAthrd + nu2 ;
  G4double xs = landa*Kc + mu + nu/Kc;

  xs = std::max(xs, 0.0);
  return xs;
}

// *********** OPT=3,4 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundNeutron::GetOpt34(G4double K)
{

  const G4double flow = 1.e-18;

  // PRECO xs for neutrons is choosen
  const G4double p0 = -312.;
  const G4double landa0 = 12.10;
  const G4double landa1=  -11.27;
  const G4double mm0 = 234.1;
  const G4double mu1 = 38.26;
  const G4double nu0 = 1.55;
  const G4double nu1 = -106.1;
  const G4double nu2 = 1280.8; 
  const G4double ra  = 0.0;

  G4double signor = 1.0;
  if(ResidualA < 40)       { signor =0.7 + ResidualA*0.0075; }
  else if(ResidualA > 210) { signor = 1. + (ResidualA-210)/250.; }

  G4double landa = landa0/ResidualAthrd + landa1;
  G4double mu = mm0*ResidualAthrd + mu1*ResidualAthrd*ResidualAthrd;
  G4double nu = nu0*ResidualAthrd*ResidualA + nu1*ResidualAthrd*ResidualAthrd + nu2;

  // JMQ very low energy behaviour corrected (problem  for A (apprx.)>60)
  if (nu < 0.) { nu = -nu; }

  G4double ec = 0.5;
  G4double ecsq = 0.25;
  G4double p = p0;
  G4double xnulam = 1.;
  G4double etest = 32.;
  //          ** etest is the energy above which the rxn cross section is
  //          ** compared with the geometrical limit and the max taken.
  //          ** xnulam here is a dummy value to be used later.

  G4double a = -2.*p*ec + landa - nu/ecsq;
  G4double b = p*ecsq + mu + 2.*nu/ec;
  G4double ecut = 0.;
  G4double cut = a*a - 4.*p*b;
  if (cut > 0.) { ecut = std::sqrt(cut); }
  ecut = (ecut-a) / (2*p);
  if (cut < 0.) { ecut -= 2.; }

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
    sig = std::max(geom,sig);
  } 
  return sig; 
}
