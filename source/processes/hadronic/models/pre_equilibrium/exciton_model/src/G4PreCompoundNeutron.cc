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

  if (OPTxs==0) { return GetOpt0( K); }
  else if( OPTxs==1 || OPTxs==2) { return GetOpt12( K); }
  else if (OPTxs==3 || OPTxs==4) { return GetOpt34( K); }
  else{
    std::ostringstream errOs;
    errOs << "BAD NEUTRON CROSS SECTION OPTION !!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;
  }
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

  G4double landa, landa0, landa1, mu, mm0, mu1,nu, nu0, nu1, nu2,xs;

  landa0 = 18.57;
  landa1 = -22.93;
  mm0 = 381.7;
  mu1 = 24.31;
  nu0 = 0.172;
  nu1 = -15.39;
  nu2 = 804.8;
  landa = landa0/ResidualAthrd + landa1;
  mu = mm0*ResidualAthrd + mu1*ResidualAthrd*ResidualAthrd;
  nu = nu0*ResidualAthrd*ResidualA + nu1*ResidualAthrd*ResidualAthrd + nu2 ;
  xs=landa*Kc + mu + nu/Kc;
  if (xs <= 0.0 ){
    std::ostringstream errOs;
    G4cout<<"WARNING:  NEGATIVE OPT=1 neutron cross section "<<G4endl;     
    errOs << "RESIDUAL: Ar=" << ResidualA << " Zr=" << ResidualZ <<G4endl;
    errOs <<"  xsec("<<Kc<<" MeV) ="<<xs <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
              }
  return xs;
}

// *********** OPT=3,4 : Kalbach's cross sections (from PRECO code)*************
G4double G4PreCompoundNeutron::GetOpt34(G4double K)
{
  G4double landa, landa0, landa1, mu, mm0, mu1,nu, nu0, nu1, nu2;
  G4double p, p0;
  G4double flow,ec,ecsq,xnulam,etest(0.),ra(0.),a,signor(1.),sig; 
  G4double b,ecut,cut,ecut2,geom,elab;

  flow = 1.e-18;

  // PRECO xs for neutrons is choosen
  p0 = -312.;
  landa0 = 12.10;
  landa1=  -11.27;
  mm0 = 234.1;
  mu1 = 38.26;
  nu0 = 1.55;
  nu1 = -106.1;
  nu2 = 1280.8; 

  if (ResidualA < 40)  { signor =0.7 + ResidualA*0.0075; }
  if (ResidualA > 210) { signor = 1. + (ResidualA-210)/250.; }
  landa = landa0/ResidualAthrd + landa1;
  mu = mm0*ResidualAthrd + mu1*ResidualAthrd*ResidualAthrd;
  nu = nu0*ResidualAthrd*ResidualA + nu1*ResidualAthrd*ResidualAthrd + nu2;

  // JMQ very low energy behaviour corrected (problem  for A (apprx.)>60)
  if (nu < 0.) { nu=-nu; }

  ec = 0.5;
  ecsq = 0.25;
  p = p0;
  xnulam = 1.;
  etest = 32.;
  //          ** etest is the energy above which the rxn cross section is
  //          ** compared with the geometrical limit and the max taken.
  //          ** xnulam here is a dummy value to be used later.

  a = -2.*p*ec + landa - nu/ecsq;
  b = p*ecsq + mu + 2.*nu/ec;
  ecut = 0.;
  cut = a*a - 4.*p*b;
  if (cut > 0.) { ecut = std::sqrt(cut); }
  ecut = (ecut-a) / (p+p);
  ecut2 = ecut;
  if (cut < 0.) { ecut2 = ecut - 2.; }
  elab = K * FragmentA / G4double(ResidualA);
  sig = 0.;
  if (elab <= ec) { //start for E<Ec 
    if (elab > ecut2) { sig = (p*elab*elab+a*elab+b) * signor; } 
  }              //end for E<Ec
  else {           //start for  E>Ec
    sig = (landa*elab+mu+nu/elab) * signor;
    geom = 0.;
    if (xnulam < flow || elab < etest) { return sig; }
    geom = std::sqrt(theA*K);
    geom = 1.23*ResidualAthrd + ra + 4.573/geom;
    geom = 31.416 * geom * geom;
    sig = std::max(geom,sig);

  }
  return sig;
}
