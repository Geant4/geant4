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
// $Id: G4eeCrossSections.cc,v 1.6 2006/06/29 19:32:42 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeCrossSections
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.10.2003
//
// Modifications:
//
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeCrossSections.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4PhysicsLinearVector.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeCrossSections::G4eeCrossSections()
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeCrossSections::~G4eeCrossSections()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeCrossSections::Initialise()
{
  MsPi = G4PionPlus::PionPlus()->GetPDGMass();
  MsPi0= G4PionZero::PionZero()->GetPDGMass();
  MsEta= 547.30*MeV;
  MsEtap=957.78*MeV;
  MsKs = G4KaonZeroLong::KaonZeroLong()->GetPDGMass();
  MsKc=G4KaonPlus::KaonPlus()->GetPDGMass();
  MsRho= 770.0*MeV;
  MsOm = 781.94*MeV;
  MsF0 = 980.0*MeV;
  MsA0 = 983.4*MeV;
  MsPhi= 1019.413*MeV;
  MsK892 = 891.66*MeV;
  MsK0892 = 896.10*MeV;
  GRho = 150.7*MeV;
  GOm = 8.41*MeV;
  GPhi = 4.43*MeV;
  GK892 = 50.8*MeV;
  GK0892 = 50.5*MeV;
  PhRho = 0.0;
  PhOm = 0.0;
  PhPhi = 155.0*degree;
  PhRhoPi = 186.0*degree;
      
  BrRhoPiG = 4.5e-4;
  BrRhoPi0G= 6.8e-4;
  BrRhoEtaG= 2.4e-4;
  BrRhoEe = 4.49e-5;
  BrOm3Pi = 0.888;
  BrOmPi0G= 0.085;
  BrOmEtaG= 6.5e-4;
  BrOm2Pi = 0.0221;
  PhOm2Pi = 90.0;
  BrOmEe = 7.07e-5;
  BrPhi2Kc = 0.491;
  BrPhiKsKl= 0.341;
  BrPhi3Pi = 0.155;
  BrPhiPi0G= 1.31e-3;
  BrPhiEtaG= 1.26e-2;
  BrPhi2Pi = 8.e-5;
  PhPhi2Pi = -20.0*degree;
  BrPhiEe = 2.99e-4;

  MsRho3 = MsRho*MsRho*MsRho;
  MsOm3  = MsOm*MsOm*MsOm;
  MsPhi3 = MsPhi*MsPhi*MsPhi;

  MeVnb = 3.8938e+11*nanobarn;
  Alpha = 1.0/137.036;

  AOmRho = 3.0;
  ARhoPRho = 0.72;
  cterm=0.;
  mssig = 600.*MeV;
  gsig = 500.*MeV;
  brsigpipi = 1.;

  msrho1450 = 1465.*MeV;
  msrho1700 = 1700.*MeV;
  grho1450 = 310.*MeV;
  grho1700 = 240.*MeV;
  arhoompi0 = 1.;
  arho1450ompi0 = 1.;
  arho1700ompi0 = 1.;
  phrhoompi0 = 0.;
  phrho1450ompi0 = pi;
  phrho1700ompi0 = 0.;
  aomrhopi0 = 1.;
  phomrhopi0 = 0.;
  arhopi0pi0g = 0.;
  aompi0pi0g = 0.;
  phrhopi0pi0g = 0.;
  phompi0pi0g = 0.;
  brrho1450ompi0 = 0.02;
  brrho1450pipi = 0.50;
  brrho1700ompi0 = 1.0;
  brrho1700pipi = 0.02;
  aphirhopi0 = 1.;
  phphirhopi0 = pi;
  arhosigg = 0.;
  phrhosigg = 0.;
  aomsigg = 0.;
  phomsigg = 0.;

  G4String w0, w1, w2;
  ph3p = 0;

  /*
  G4double emin, emax;
  G4int nbins;
  const G4String fname = "wrhopi.wid"; 
  ifstream  fi(fname.c_str()); 
  fi >> w0 >> nbins >> w1 >> emin >> w2 >> emax;
  emin *= MeV;
  emax *= MeV;
  ph3p = new G4PhysicsLinearVector(emin,emax,nbins);
  G4int nlines = nbins/5;
  G4double s0, s1, s2, s3, s4;
  for(G4int i=0; i<nlines; i++) {
    fi >> s0 >> s1 >> s2 >> s3 >> s4;
    ph3p->PutValue(5*i, s0);
    ph3p->PutValue(5*i + 1, s1);
    ph3p->PutValue(5*i + 2, s2);
    ph3p->PutValue(5*i + 3, s3);
    ph3p->PutValue(5*i + 4, s4);
  }
  fi.close();
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::CrossSection2pi(G4double e) 
{
 
  complex<G4double> xr(cos(PhRho),sin(PhRho));
  complex<G4double> xo(cos(PhOm2Pi),sin(PhOm2Pi));
  complex<G4double> xf(cos(PhPhi2Pi),sin(PhPhi2Pi));

  G4double s = e*e;
  complex<G4double> drho = DpRho(e);
  complex<G4double> dom  = DpOm(e);
  complex<G4double> dphi = DpPhi(e);

  complex<G4double> amp = 
      sqrt(Width2p(s,MsRho,GRho,1.0,MsPi)*MsRho3*BrRhoEe*GRho)*xr/drho
    + sqrt(Width2p(s,MsOm,GOm,BrOm2Pi,MsPi)*MsOm3*BrOmEe*GOm)*xo/dom
    + sqrt(Width2p(s,MsPhi,GPhi,BrPhi2Pi,MsPi)*MsPhi3*BrPhiEe*GPhi)*xf/dphi;

  G4double cross = 12.0*pi*MeVnb*norm(amp)/(e*s);

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::Width2p(G4double s, G4double mres, 
                                    G4double g, G4double br, G4double m) 
{
  G4double m2    = 4.0*m*m;
  G4double s0    = mres*mres;
  G4double f     = (s - m2)/(s0 - m2);
  if(f < 0.0)  f = 0.0;
  return g*br*sqrt(f)*f*s0/s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::Width3p(G4double s, G4double mres, 
                                    G4double g, G4double br) 
{
  G4double w = PhaseSpace3p(sqrt(s));
  G4double w0= PhaseSpace3p(mres);
  G4double x = g*br*w/w0;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::PhaseSpace3p(G4double e) 
{
  
  //  G4bool b;
  //  G4double x = ph3p->GetValue(e, b);
  G4double x = 1.0; 
  G4double emev = e/MeV;
  G4double y = 414.12/emev;
  x *= pow(e/MsOm, 5.0) * pow(emev*0.1, 3.0)*(1.0 - y*y);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::WidthPg(G4double s, G4double mres, 
                                    G4double g, G4double br, G4double m)
{
  G4double m2    = m*m;
  G4double s0    = mres*mres;
  G4double f     = (s - m2)*mres/((s0 - m2)*sqrt(s));
  if(f < 0.0)  f = 0.0;
  return g*br*f*f*f;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::WidthRho(G4double e) 
{
  G4double w = Width2p(e*e, MsRho, GRho, 1.0, MsPi);
  return w;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::WidthOm(G4double e) 
{
  G4double s = e*e;
  G4double w = (Width3p(s, MsOm, GOm, BrOm3Pi) +
                WidthPg(s, MsOm, GOm, BrOmPi0G, MsPi0) +
                WidthPg(s, MsOm, GOm, BrOmEtaG, MsEta) +
                Width2p(s, MsOm, GOm, BrOm2Pi, MsPi)) /
       (BrOm3Pi+BrOmPi0G+BrOmEtaG+BrOm2Pi);
  return w;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeCrossSections::WidthPhi(G4double e) 
{
  G4double s = e*e;
  G4double w = (Width3p(s, MsPhi, GPhi, BrPhi3Pi) +
                WidthPg(s, MsPhi, GPhi, BrPhiPi0G, MsPi0) +
                WidthPg(s, MsPhi, GPhi, BrPhiEtaG, MsEta) +
                Width2p(s, MsPhi, GPhi, BrPhi2Kc, MsKc) +
                Width2p(s, MsPhi, GPhi, BrPhiKsKl, MsKs)) /
        (BrPhi3Pi+BrPhiPi0G+BrPhiEtaG+BrPhi2Kc+BrPhiKsKl);
  return w;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

complex<G4double> G4eeCrossSections::DpRho(G4double e) 
{
  complex<G4double> d(MsRho*MsRho - e*e, -e*WidthRho(e));
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

complex<G4double> G4eeCrossSections::DpOm(G4double e) 
{
  complex<G4double> d(MsOm*MsOm - e*e, -e*WidthOm(e));
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

complex<G4double> G4eeCrossSections::DpPhi(G4double e) 
{
  complex<G4double> d(MsPhi*MsPhi - e*e, -e*WidthPhi(e));
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
