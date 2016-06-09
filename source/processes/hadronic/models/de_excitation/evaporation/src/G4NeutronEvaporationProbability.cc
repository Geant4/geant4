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
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03-09-2008 J.M. Quesada for external choice of inverse cross section option
// 17-11-2010 V.Ivanchenko integer Z and A

#include "G4NeutronEvaporationProbability.hh"
#include "G4SystemOfUnits.hh"

G4NeutronEvaporationProbability::G4NeutronEvaporationProbability() :
    G4EvaporationProbability(1,0,2,&theCoulombBarrier) // A,Z,Gamma,&theCoulombBarrier
{
  ResidualA = ResidualZ = theA = theZ = FragmentA = 0;
  ResidualAthrd = FragmentAthrd = 0.0;
}

G4NeutronEvaporationProbability::~G4NeutronEvaporationProbability()
{}

G4double G4NeutronEvaporationProbability::CalcAlphaParam(const G4Fragment & fragment) 
{ 
  return 0.76+2.2/fG4pow->Z13(fragment.GetA_asInt()-GetA());
}
	
G4double G4NeutronEvaporationProbability::CalcBetaParam(const G4Fragment &  fragment) 
{ 
  return (2.12/fG4pow->Z23(fragment.GetA_asInt()-GetA()) - 0.05)*MeV/
    CalcAlphaParam(fragment); 
}

////////////////////////////////////////////////////////////////////////////////////
//J. M. Quesada (Dec 2007-June 2008): New inverse reaction cross sections 
//OPT=0 Dostrovski's parameterization
//OPT=1,2 Chatterjee's paramaterization 
//OPT=3,4 Kalbach's parameterization 
// 
G4double 
G4NeutronEvaporationProbability::CrossSection(const  G4Fragment & fragment, G4double K)
{ 
  theA=GetA();
  theZ=GetZ();
  ResidualA=fragment.GetA_asInt()-theA;
  ResidualZ=fragment.GetZ_asInt()-theZ; 
  
  ResidualAthrd=fG4pow->Z13(ResidualA);
  FragmentA=fragment.GetA_asInt();
  FragmentAthrd=fG4pow->Z13(FragmentA);

  if (OPTxs==0) {std::ostringstream errOs;
    errOs << "We should'n be here (OPT =0) at evaporation cross section calculation (neutrons)!!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;}
  else if( OPTxs==1 ||OPTxs==2) return GetOpt12( K);
  else if (OPTxs==3 || OPTxs==4)  return GetOpt34( K);
  else{
    std::ostringstream errOs;
    errOs << "BAD NEUTRON CROSS SECTION OPTION AT EVAPORATION!!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    return 0.;
  }
}
 
//********************* OPT=1,2 : Chatterjee's cross section ***************
//(fitting to cross section from Bechetti & Greenles OM potential)

G4double G4NeutronEvaporationProbability::GetOpt12(G4double K)
{
  G4double Kc=K;

  // Pramana (Bechetti & Greenles) for neutrons is chosen 

  // JMQ  xsec is set constat above limit of validity
  if (K > 50*MeV) { Kc = 50*MeV; }

  G4double landa, landa0, landa1, mu, mum0, mu1,nu, nu0, nu1, nu2,xs;

  landa0 = 18.57;
  landa1 = -22.93;
  mum0 = 381.7;
  mu1 = 24.31;
  nu0 = 0.172;
  nu1 = -15.39;
  nu2 = 804.8;
  landa = landa0/ResidualAthrd + landa1;
  mu = mum0*ResidualAthrd + mu1*ResidualAthrd*ResidualAthrd;
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
G4double G4NeutronEvaporationProbability::GetOpt34(G4double K)
{
  G4double landa, landa0, landa1, mu, mum0, mu1,nu, nu0, nu1, nu2;
  G4double p, p0;
  G4double flow,ec,ecsq,xnulam,etest(0.),ra(0.),a,signor(1.),sig; 
  G4double b,ecut,cut,ecut2,geom,elab;

  //safety initialization
  landa0=0;
  landa1=0;
  mum0=0.;
  mu1=0.;
  nu0=0.;
  nu1=0.;
  nu2=0.;
  p0=0.;

  flow = 1.e-18; 

  // PRECO xs for neutrons is choosen
  p0 = -312.;
  landa0 = 12.10;
  landa1=  -11.27;
  mum0 = 234.1;
  mu1 = 38.26;
  nu0 = 1.55;
  nu1 = -106.1;
  nu2 = 1280.8; 

  if (ResidualA < 40)  { signor =0.7 + ResidualA*0.0075; }
  if (ResidualA > 210) { signor = 1. + (ResidualA-210)/250.; }
  landa = landa0/ResidualAthrd + landa1;
  mu = mum0*ResidualAthrd + mu1*ResidualAthrd*ResidualAthrd;
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

