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
//
// G4 Tools program: Glauber Hadron-Nuclear Cross Sections
// .....................................................
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Jan-2005
// 
//=====================================================================

//#define debug
//#define bestest
//#define integrcs
#define eldiffcs

#include "globals.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include "G4ios.hh"
#include "G4QBesIKJY.hh"
#include "G4QPDGCode.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4NucleiProperties.hh"
//#include <CLHEP/GenericFunctions/Bessel.hh>

// Integrated cross-sections
G4double TotalCrSec, InelCrSec, ProdCrSec, ElasticCrSec, QuasyElasticCrSec; 
// Parameters of hadron-nucleon interaction
G4double HadrTot, HadrSlope, HadrReIm, DDSect2, DDSect3, Tot00;
// CHIPS parameters of hadron-proton & hadron-neutron interaction (Re/Im is common)
G4double HadrTElN, HadrTElP, HadrSlpN, HadrSlpP, HadrPexN, HadrPexP;

// CHIPS Parameters (only for protons, nn is necessary for neutrons np+nn)
void CHIPSParameters(G4int iPDG, G4double HadrMoment) // HadrMoment is in MeV
{
  G4double p=HadrMoment/GeV;                   // Momentum in GeV
  G4double p2=p*p;                             // Squared momentum in GeV^2
  G4double sp=std::sqrt(p);
  G4double p2s=p2*sp;
  G4double p3=p2*p;
  G4double p4=p2*p2;
  G4double p5=p4*p;
  G4double lp=std::log(p);
  G4double d2=(lp-5.)*(lp-5.);
  G4double p8=p4*p4;
  if (iPDG==2212)
  {
    HadrTElN=12./(p2s+.05*p+.0001/sp)+.35/p+(6.75+.14*d2+19./p)/(1.+.6/p4); // pn (mb)
    HadrTElP=2.91/(p2s+.0024)+(6.75+.14*d2+23./p)/(1.+.4/p4);               // pp (mb)
    HadrSlpN=(7.2+4.32/(p8+.012*p3))/(1.+2.5/p4);                           // pn (GeV^-2) 
    HadrSlpP=8.*std::pow(p,.055)/(1.+3.64/p4);                              // pp (GeV^-2) 
    HadrPexN=(6.75+.14*d2+13./p)/(1.+.14/p4)+.6/(p4+.00013);                // pn (mb/sr)
    HadrPexP=(74.+3.*d2)/(1.+3.4/p5)+(.2/p2+17.*p)/(p4+.001*sp);            // pn (mb/sr)
    HadrReIm=-.54+.0954*lp+.975/(p2+.09865/p);           // no unit (the same for pp & pn)
    //G4cout<<"CHIPSPar:p="<<p<<",N="<<HadrTotN<<",P="<<HadrTotN<<",RI="<<HadrReIm<<G4endl;
  }     
  else G4cout<<"CHIPSPar: Only proton projectile (2212) is implemented PDG="<<iPDG<<G4endl;
}

// Dubna Parameters
G4int CalculateParameters(G4int iPDG, G4double HadrMoment) // HadrMoment is in MeV
{
  static const G4double mN=.938;   // Dubna's mass of the nucleon in nuclei in GeV
  //static const G4double mN=.9315;  // AtomicMassUnit (mass of the nucleon in nuclei) in GeV
  static const G4double dmN=mN+mN; // Doubled mass of the nucleon in nuclei in GeV
  static const G4double mN2=mN*mN; // Squared mass of the nucleon in nuclei in GeV^2
  G4int iHadron=-1;
  if(iPDG==2212||iPDG==2112||iPDG==3122||iPDG==3222||iPDG==3112||
     iPDG==3212||iPDG==3312||iPDG==3322||iPDG==3334) iHadron = 0;
  else if(iPDG==-2212||iPDG==-2112||iPDG==-3122||iPDG==-3222||
          iPDG==-3112||iPDG==-3212||iPDG==-3312||iPDG==-3322||
          iPDG==-3334) iHadron = 1;       // Anti-nucleons, Anti-hyperons
  else if(iPDG== 211)  iHadron = 2;       // Pi+
  else if(iPDG==-211)  iHadron = 3;       // Pi-
  else if(iPDG== 321)  iHadron = 4;       // K+
  else if(iPDG==-321)  iHadron = 5;       // K- @@ What about K0?
  else
  {
    G4cout<<"CalculateParameters: PDG= "<<iPDG<<" is not supported"<<G4endl;
    return iHadron;
  }
  G4double mHadr      = G4QPDGCode(iPDG).GetMass()/1000.; // Hadron Mass in GeV
  G4double mHadr2     = mHadr*mHadr;               // Squared Hadron Mass in GeV
           HadrMoment/= 1000.;                     // Momentum in GeV (input parameter)
  G4double HadrMoment2= HadrMoment*HadrMoment;     // Squared Momentum in GeV^2
  G4double HadrEnergy2= mHadr2+HadrMoment2;        // Tot energy in GeV
  G4double HadrEnergy = std::sqrt(HadrEnergy2);    // Tot energy in GeV (global parameter)
  G4double sHadr      = HadrEnergy*dmN+mN2+mHadr2; // s in GeV^2
  G4double sqrS       = std::sqrt(sHadr);          // W in GeV
  //G4cout<<"GHAD: E="<<HadrEnergy<<",dM="<<dmN<<",sN="<<mN2<<",sH="<<mHadr2<<G4endl;
  switch (iHadron)
  {
  case 0:                                      // =====>>> proton or hyperons (the same?)
    //G4double Delta=1;                                                   // LowEnergy corr
    //if(HadrEnergy<40.) Delta = 0.916+0.0021*HadrEnergy;
    HadrTot = 5.2+5.2*std::log(HadrEnergy)+51*std::pow(HadrEnergy,-0.35); // mb
    HadrSlope = 6.44+0.88*std::log(sHadr)-1;                              // GeV^-2 
    HadrReIm  = 0.13*std::log(sHadr/350)*std::pow(sHadr,-0.18);           // no unit
    //G4cout<<"GHAD:S="<<sHadr<<",T="<<HadrTot<<",R="<<HadrReIm<<",B="<<HadrSlope<<G4endl;
    DDSect2   = 11.;                                                      // mb*GeV^-2
    DDSect3   = 3.;                                                       // mb*GeV^-2
    // Hyperon corrections
    if(iPDG==3122||iPDG==3222||iPDG==3112||iPDG==3212) //Lambda,Sig0+-
    {
      HadrTot   *= 0.80;
      HadrSlope *= 0.85;
    }
    else if(iPDG==3312||iPDG==3322)            // ---> Xi0-
    {
      HadrTot   *= 0.70;
      HadrSlope *= 0.75;
    }           
    else if(iPDG==3334)                        // ---> Omega-
    {
      HadrTot   *= 0.60;
      HadrSlope *= 0.65;
    }
    break;
  case 1:                                      // =====>>> antiproton and antihyperons
    HadrTot = 5.2+5.2*std::log(HadrEnergy)+123.2*std::pow(HadrEnergy,-0.5); //  mb
    HadrSlope = 8.32+0.57*std::log(sHadr);     // GeV^-2 
    if(HadrEnergy<1000.) HadrReIm = 0.06*(sqrS-2.236)*(sqrS-14.14)*std::pow(sHadr,-0.8);
    else HadrReIm  = 0.6*std::log(sHadr/350)*std::pow(sHadr,-0.25);
    DDSect2 = 11.;                             //mb*GeV-2 @@ the same as in case 0
    DDSect3 = 3.;                              //mb*GeV-2 @@ the same as in case 0
    // Hyperon corrections
    if(iPDG==-3122||iPDG==-3222||iPDG==-3112||iPDG==-3212)//AntiLam,Sig
    {
      HadrTot   *=0.75;
      HadrSlope *=0.85;
    }
    else if(iPDG==-3312||iPDG==-3322)          // =====>>> AntiXi0-
    {
      HadrTot   *=0.65;
      HadrSlope *=0.75;
    }
    if(iPDG==-3334)                            // ======>>>AntiOmega-
    {
      HadrTot   *=0.55;
      HadrSlope *=0.65;
    }
    break;
  case 2:                                      // =====>>>  pi plus
    if(HadrMoment>2.) HadrTot = 10.6+2.*std::log(HadrEnergy)+25*std::pow(HadrEnergy,-0.43);// mb
    else HadrTot = 40-50*(HadrMoment-1.5)*(HadrMoment-1.7);
    HadrSlope = 7.28+0.245*std::log(sHadr);                                // GeV^-2 
    HadrReIm  = 0.2*std::log(sHadr/100)*std::pow(sHadr,-0.15);             // no dim
    DDSect2   = 4.6;                                                       // mb*GeV^-2
    DDSect3   = 1.33;                                                      // mb*GeV^-2
    break;
  case 3:                                      // =====>>>  pi minus
    HadrTot = 10.6+2*std::log(HadrEnergy)+30*std::pow(HadrEnergy,-0.43);   // mb @@ E->inf?
    if(HadrMoment<1.399) HadrTot = HadrTot+21.0/0.4*(1.4-HadrMoment);      // A huge jump ?
    HadrSlope = 7.28+0.245*std::log(sHadr);                                // GeV-2 
    HadrReIm  = 0.2*std::log(sHadr/100)*std::pow(sHadr,-0.15);
    DDSect2   = 4.6;                           //mb*GeV-2 @@ the same as in case 2
    DDSect3   = 1.33;                          //mb*GeV-2
    break;
  case 4:                                      // =====>>> K plus
    HadrTot = 10.6+1.8*std::log(HadrEnergy)+9.*std::pow(HadrEnergy,-0.55); // mb 
    if(HadrEnergy>100) HadrSlope = 15.;              // Jump at 100 GeV (?)
    else HadrSlope = 1.0+1.76*std::log(sHadr)-2.84*std::pow(sHadr,-0.5);     // GeV-2
    //else HadrSlope = 5.28+1.76*std::log(sHadr)-2.84*std::pow(sHadr,-0.5);    // GeV-2
    HadrReIm = 0.4*(sHadr-20)*(sHadr-150)*std::pow(sHadr+50,-2.1);
    DDSect2   = 3.5;                           //mb*GeV-2
    DDSect3   = 1.03;                          //mb*GeV-2
    break;
  case 5:                                      // ======>>> K minus
    HadrTot = 10+1.8*std::log(HadrEnergy)+25*std::pow(HadrEnergy,-0.5);      // mb 
    HadrSlope = 6.98+0.127*std::log(sHadr);           // GeV-2 
    //if(HadrEnergy<8) HadrReIm = 0.7;
    HadrReIm  = 0.4*(sHadr-20)*(sHadr-20)*std::pow(sHadr+50,-2.1);
    DDSect2   = 3.5;                           //mb*GeV-2
    DDSect3   = 1.03;                          //mb*GeV-2
    break;
  }     
  return iHadron;
}

void  CalculateIntegralCS(G4int  Anuc, G4double HadrEnergy)
{
  static const G4double eps = 0.0001;          // Accuracy of calculations
  static const G4double mb2G2 = 2.568;         // Transformation from mb to GeV^-2
  static const G4double incrf = 10;            // Increase factor of radius for intergation
  static const G4double m2G10 = incrf*mb2G2;   // Big conversion factor
  static const G4double dR0   = 1.06+1.06;     // Doubled nominal R0
  G4double An      = Anuc;                     // A float
  G4double Stot    = HadrTot*mb2G2;            // Total hadron-nucleon CS in GeV-2
  G4double Bhad    = HadrSlope;                // Diffraction cone nH In GeV-2
  G4double Asq     = 1+HadrReIm*HadrReIm;      // |M|^2/(ImM)^2

  G4double R0      = std::sqrt(0.99);          // in fermi (strange value? 1.07?)

  if (Anuc == 58) R0 = std::sqrt(0.6);         // Personal correction for Fe
  if (Anuc >20)   R0 = std::sqrt((35.34+0.5*Anuc)/(40.97+Anuc));
  if (Anuc == 16) R0 = std::sqrt(0.75);        // Personal correction for O16
  if (Anuc >10)   R0 = std::sqrt(0.84);        // The same R0 for A=11-20 (except O16)
  if (Anuc ==  4) R0 *= 1.2;                   // Personal correction for Helium

  G4double Rnucl   = R0*std::pow(static_cast<double>(Anuc),.3333333); // in fermi
  G4double Rnuc2   = Rnucl*Rnucl*m2G10;        // in GeV^-2 (10 -> to make it >> R)
  G4double RB      = Rnuc2+Bhad;               // Effective increase of R2 (Convolution?)
  G4double R2B     = RB+Bhad;                  // Over-convolution (?)
  G4double Delta   = Stot/R2B/twopi;           // Integration step for Total & Inelastic CS
  G4double Delt2   = Delta+Delta;              // Two delta
  G4double Delt3   = Stot*Asq*R2B/RB/Bhad/16/pi; // Integration step for Production CS

  G4double Tot0    = 0.;                       // =SUM[(-Delta)^(i-1)*(C^A_i)/i]
  G4double Inel0   = 0.;                       // =SUM[(-Delta)^(i-1)*(2A)!/(2A-i)!/i!/i]
  G4double N       = -1./Delta; 
  G4double N1      = N;
  G4double N3      = -1./Delt2;
  G4double Prod0   = 0.;
  G4int    Ai      = Anuc+1;                   // Working value to avoid multiplication
  G4double RBi     = 0.;                       // Working value to avoid multiplication
  for (G4int i=1; i<= Anuc; i++)               // ==> Loop over nucleons
  {
    Ai--;                                      // Working value to avoid multiplication
    RBi   += RB;                               // Working value to avoid multiplication
    G4double Di=Delta/i;                       // Working value to avoid repetition
    N     *= -Di*Ai;                           // (-Delta)^(i-1)*A!/(A-i)!/i! (C^A_i)
    N1    *= -Di*(Anuc+Ai);                    // (-Delta)^(i-1)*(2A)!/(2A-i)!/i!
    N3    *= -Delt2*Ai/i;                      // (-2*Delt)^(i-1)*A!/(A-i)!/i! (C^A_i)
    Tot0  += N/i;                              // =SUM[(-Delta)^(i-1)*(C^A_i)/i]
    Inel0 += N1/i;                             // =SUM[(-Delta)^(i-1)*(2A)!/(2A-i)!/i!/i]
    //G4double N2    = 1.;
    G4double N4    = 1.;
    //G4double Inel1 = 0.;
    //G4double Inel2 = 0.;
    G4double Prod1 = 0.;
    G4double Bhl   = 0.;
    for (G4int l=0; l<= i; l++)                // ==> Loop over nucleon-nucleon pairs (0?)
    {
      //Inel1 += N2/(i+l);                 // SUM[(-Delta)^(l-1)*i!(i-1)!/(l-1)!/l!/(i+l)!]
      Prod1 += N4*RB/(RBi+Bhl);                //
      //N2    *= -Delt*(i-l)/(l+1);              // (-Delta)^l*i!/l!/(l+1)!
      N4    *= -Delt3*(i-l)/(l+1);
      Bhl   += Bhad;
    }
    //Inel2 += Inel1*N3;
    Prod0   += Prod1*N3;
    if(std::fabs(N1/i/Inel0) < eps)  break;
  }

  Tot0             *= HadrTot;                 // Multiply by the hN cross section
  Inel0            *= HadrTot/2;               // Multiply by a half of hN cross section
  //Inel2            *= HadrTot;                 // Multiply by the hN cross section
  Prod0            *= HadrTot;                 // Multiply by the hN cross section
  Tot00             = Tot0;                    // Remember because it will be changed
  G4double ak       = Rnuc2*twopi/Stot;        // cross-section ratio for integration
  G4double Ank      = An/ak;                   // Working value to avoid repetition
  G4double rak      = 1./ak;                   // Working value to avoid repetition
  G4double ak1      = 1.-rak/4;                // Working value to avoid repetition
  G4double Ank1     = Ank*ak1;                 // Working value to avoid repetition
  G4double DDSect1  = (DDSect2+DDSect3*std::log(dR0*HadrEnergy/Rnucl/std::sqrt(m2G10)/4));
  G4double Dtot     = 8*pi*ak/HadrTot*(1.-(1.+Ank)*std::exp(-Ank))*DDSect1/mb2G2;//TotCSCor
  G4double bk       = (1.-rak)/Stot/ak1;       // Working value to avoid power(,2)
  G4double bd       = bk*bk*DDSect1*(1.-(1.+Ank1)*std::exp(-Ank1))*Rnuc2; // Working
  G4double Dprod    = bd*4*pi2*mb2G2;          // Correction for the Production & Inel CS
  TotalCrSec        = Tot0-Dtot;               // Total Cross Section         (primary)
  InelCrSec         = Inel0-Dprod;             // Inelastic Cross Section     (primary)
  //InelCrSec1        = Inel2;                 // Direct calculation without correction
  ProdCrSec         = Prod0-Dprod;             // Production Cross Section    (primary)
  ElasticCrSec      = TotalCrSec-InelCrSec;    // Elastic Cross Section       (secondary)
  QuasyElasticCrSec = InelCrSec-ProdCrSec;     // Quasy-Elastic Cross Section (secondary)
}

#ifdef eldiffcs
// Functions for calculation of differential elastic cross-sections
// Static globals
static const G4double  rAmax = 2.5;      // Factor of maximum radius
static const G4int     NptB  = 500;      // A#of steps in the impact parameter integration
static const G4int     nFact = 253;      // Maximum A (for Factorials and Factors)

// Globals
G4double r0 = 1.1;                       // Wood-Saxon's factor
G4double Mnoj[nFact], Factorials[nFact]; // Factors and factorials
G4double R1, R2, Pnucl, Aeff;            // Nuclear Parameters
G4double InCoh;                          // Incoherent quasi-elastic cross section

G4double binom(G4int N, G4int M)
{
  if(N>170 && N>M && M)
		{
    G4double fN=N;
    G4double fM=M;
    G4double fD=N-M;
    G4double man=N*std::log(fN)-M*std::log(fM)-fD*std::log(fD)-1.;
    //G4cout<<"G4GCS::binom: N="<<N<<", M="<<M<<", C="<<man<<", E="<<std::exp(man)<<G4endl;
    return std::exp(man);
		}
  else if(N>1 && N>M && M) return Factorials[N]/Factorials[M]/Factorials[N-M];
  return 1.;
}

// Initialize factorials and combinatoric coefficients
void Initialize()
{
  //static const G4double  sat = 3.03;      // Stirling saturation of Dubna
  G4double Val = Factorials[0] = Mnoj[0] = 1.;
  for (G4int i=1; i<nFact; i++)
  {
    G4double Si = 0.;
             Val *= i;                    // Calculation of the factorial values
    Factorials[i] = Val;
    for (G4int l = 0; l<=i; l++)
    {
      G4double Sm = 1.;                   // sum of reversed binom coefficients (?)
      for (G4int m=1; m<=i-l; m++) Sm += 1./binom(i-l,m);
      //if(i>169)G4cout<<"G4GCS::Init:S="<<Sm<<", b("<<i<<","<<l<<")="<<binom(i,l)<<G4endl;
      Si += Sm/binom(i,l);
    } //  End of l LOOP
    Mnoj[i] = Si;
    //G4double d=Si-sat;
    //G4double r=d/sat;
    //if(i>100)G4cout<<"G4GCS::Init:i="<<i<<":"<<Si<<"-"<<sat<<" = "<<d<<", r="<<r<<G4endl;
  } // End of i LOOP
} // End of Initialization

// Dubna paranmeters of nuclei
void CalculateNuclearParameters(G4int Anuc)  // R1,R2,Pnuc,Aeff
{//  ==============================================
  G4double A=Anuc;
  // Personal corrections fore some nuclei
  if(Anuc == 4)                                      // Personal correction for He (!)
  {
    R1    = 5.5;
    R2    = 3.7;
    Pnucl = 0.4;   
    Aeff  = 0.87;
  }
  else if(Anuc == 9)                                 // Personal correction for Be (!)
  {
    R1    = 9.0;
    R2    = 7.0;
    Pnucl = 0.190;
    Aeff  = 0.9;
  }
  //else if(Anuc == 11)                              // Personal correction for B (!)
  //{
  //  R1    = 10.8;
  //  R2    = 7.5;
  //  Pnucl = 0.85;
  //  Aeff  = 1.2;
  //}
		//else if(Anuc == 12)                              // Personal correction for C (!)
  //{
  //   R1    = 9.336;
  //   R2    = 5.63;
  //   Pnucl = 0.197;
  //   Aeff  = 1.0;
  //}
  else if(Anuc == 16)                                // Personal correction for O (!)
  {
    R1    = 10.50;
    R2    = 5.5;
    Pnucl = 0.7;
    Aeff  = 0.98;
    //R1    = 11.3;
    //R2    = 2.5;
    //Pnucl = 0.75;
    //Aeff  = 0.9;
  }
  else if(Anuc == 58)                                // Personal correction for Ni (!)
  {
    R1    = 15.0;
    R2    = 9.9;
    Pnucl = 0.45;
    Aeff  = 0.85;
  }
  else if(Anuc == 90)                                // Personal correction for Zr (!)
  {
    //R1    = 16.5;
    //R2    = 11.62;
    //Pnucl = 0.4;
    //Aeff  = 0.9;
    R1    = 16.5;
    R2    = 11.62;
    Pnucl = 0.4;
    Aeff  = 0.7;
  }
  else if(Anuc == 208)                               // Personal correction for Pb (!)
  {  
    //      R1 = 20.73; R2 = 15.74.
    //R1    = 4.1408*std::pow(A,0.3018);
    //R2    = 3.806*std::pow(A-10.068,0.2685);
    //Pnucl = 0.9;
    //Aeff  = 1.1;
    R1    = 19.5;
    R2    = 15.74;
    Pnucl = 0.4;
    Aeff  = 0.7;
  }
  else
  {
    R1 = 4.45*std::pow(A-1.,.309);                  // First diffractional radius
    if(Anuc == 28) R1 *= .955;                       // Personal correction for Si (?!)
    R2    = 2.3*std::pow(A,.36);                    // Second diffraction radius
    Pnucl = 0.176+(.00167+.00000869*Anuc)*Anuc;      // Momentum of mean Fermi motion
    Aeff  = 0.9;                                     // Effective screaning
  }
}

// Independent transition method calculation. Has problems for A>95 (short cut).
G4double DiffElasticCS(G4int hPDG, G4double HadrMom, G4int A, G4double aQ2) // All in MeV
{//      =================================================================
  static const G4double eps = 0.000001;                     // Accuracy of calculations
  static const G4double mb2G2 = 2.568;                      // Transform from mb to GeV^-2
  static const G4double piMG  = pi/mb2G2;                   // PiTrans from GeV^-2 to mb
  static const G4double dR0   = 1.06+1.06;                  // Doubled nominal R0
  G4double hMom       = HadrMom/1000.;                      // Momentum (GeV, inputParam)
		G4double MassH      = G4QPDGCode(hPDG).GetMass()/1000.;   // Hadron Mass in GeV
  G4double MassH2     = MassH*MassH;
  G4double HadrEnergy = std::sqrt(hMom*hMom+MassH2);        // Tot energy in GeV
		//G4cout<<"G4GCS::DiffElasticCS: PGG_proj="<<hPDG<<", A_nuc= "<<A<<G4endl;
  if(A==2 || A==3) G4Exception("G4GCS: This model does not work for nuclei with A=2, A=3");
  if(A>252) G4Exception("G4GCS: This nucleus is too heavy for the model !!!");
  //if(HadrEnergy-MassH<1.) G4Exception("Kin energy is too small for the model (T>1 GeV)");
  CalculateParameters(hPDG, HadrMom);                   
  CalculateIntegralCS(A, HadrEnergy);
  CalculateNuclearParameters(A);
  G4double Q2 = aQ2/1000/1000;                              // in GeV^2
		//MassN           = A*0.938;                              // @@ This is bad!
		G4double MassN    = A*0.9315;                             // @@ Atomic Unit is bad too
  G4double MassN2   = MassN*MassN;
  G4double S        = (MassN+MassN)*HadrEnergy+MassN2+MassH2;// Mondelststam s
  G4double EcmH     = (S-MassN2+MassH2)/2/std::sqrt(S);     // CM energy of a hadron
  G4double CMMom    = std::sqrt(EcmH*EcmH-MassH2);          // CM momentum

  G4double Stot     = HadrTot*mb2G2;                        // in GeV^-2
  G4double Bhad     = HadrSlope;                            // in GeV^-2
  G4double Asq      = 1+HadrReIm*HadrReIm;                  // |M|^2/(ImM)^2
  G4double Rho2     = std::sqrt(Asq);                       // M/ImM
  G4double R12      = R1*R1;
  G4double R22      = R2*R2;
  G4double R12B     = R12+Bhad+Bhad;
  G4double R22B     = R22+Bhad+Bhad;
  G4double R12Ap    = R12+20;
  G4double R22Ap    = R22+20;
  G4double R13Ap    = R12*R1/R12Ap;
  G4double R23Ap    = R22*R2/R22Ap*Pnucl;
  G4double R23dR13  = R23Ap/R13Ap;
  G4double R12Apd   = 2./R12Ap;
  G4double R22Apd   = 2./R22Ap;
  G4double R122Apd  = (R12Apd+R22Apd)/2;
  G4double dR0H     = dR0*HadrEnergy/4;
  G4double DDSec1p  = DDSect2+DDSect3*std::log(dR0H/R1);
  G4double DDSec2p  = DDSect2+DDSect3*std::log(dR0H/std::sqrt((R12+R22)/2));
  G4double DDSec3p  = DDSect2+DDSect3*std::log(dR0H/R2);
  G4double Norm     = (R12*R1-Pnucl*R22*R2)*Aeff;    // Some questionable norming
  G4double R13      = R12*R1/R12B;
  G4double R23      = Pnucl*R22*R2/R22B;
  G4double norFac   = Stot/twopi/Norm;               // totCS (in GeV^-2) is used
  G4double Unucl    = norFac*R13;
  G4double UnuclScr = norFac*R13Ap;
  G4double SinFi    = HadrReIm/Rho2;
  G4double FiH      = std::asin(SinFi);
  G4double N        = -1.;
  G4double N2       = R23/R13;
  G4double ImElA0   = 0.;
  G4double ReElA0   = 0.;
  G4double Tot1     = 0.;
  for(G4int i=1; i<=A; i++)                          // @@ Make separately for n and p
  {
    N              *= -Unucl*(A-i+1)*Rho2/i;         // Includes total cross-section
    G4double N4     = 1.;
    G4double Prod1  = std::exp(-Q2*R12B/i/4)*R12B/i; // Includes the slope
    G4double medTot = R12B/i;
    for(G4int l=1; l<=i; l++)
    {
      G4double exp1 = l/R22B+(i-l)/R12B;
      N4           *= -(i-l+1)*N2/l;
      G4double Nexp = N4/exp1;
      Prod1        += Nexp*std::exp(-Q2/exp1/4);
      medTot       += Nexp;
    }
    G4double FiHi   = FiH*i;
    G4double nCos   = N*std::cos(FiHi);
    ReElA0 += Prod1*N*std::sin(FiHi);
    ImElA0 += Prod1*nCos;
    Tot1           += medTot*nCos;
    if(std::abs(Prod1*N/ImElA0) < eps) break;
  }
  ImElA0   *= piMG;                                  // The amplitude in mb
  ReElA0   *= piMG;                                  // The amplitude in mb
  Tot1             *= piMG+piMG;
  G4double N1p      = 1.;
  G4double R13Ap1   = DDSec1p*R13Ap*R13Ap/2;
  G4double R22Ap1   = DDSec2p*R13Ap*R23Ap;
  G4double R23Ap1   = DDSec3p*R23Ap*R23Ap*R22Ap/2;
  G4double R13Ap2   = R13Ap1*R12Ap/4;
  G4double R22Ap2   = R22Ap1/R122Apd/2;
  G4double R23Ap2   = R23Ap1*R22Ap/4;
  G4double Q28      = Q2/8;
  G4double Q24      = Q28+Q28;
  G4double Din1     = R13Ap2*std::exp(-Q28*R12Ap)-R22Ap2*std::exp(-(Q28+Q28)/R122Apd)+
				                  R23Ap2*std::exp(-Q28*R22Ap);   // i=0 start value
  G4double DTot1    = R13Ap2-R22Ap2+R23Ap2;          // i=0 start value
  if (A<96)
		{
    for(G4int i=1; i<A-1; i++)
    {
      N1p              *= -UnuclScr*(A-i-1)/i*Rho2;
      G4double N2p      = 1.;
      G4double Din2     = 0.;
      G4double DmedTot  = 0.;
      G4double BinCoeff = 1.;                        // Start for Binomial Coefficient
      for(G4int l = 0; l<=i; l++) 
      {
        if(l) BinCoeff *= (i-l+1)/l;
        G4double exp1   = l/R22B+(i-l)/R12B;
        G4double exp1p  = exp1+R12Apd;
        G4double exp2p  = exp1+R122Apd;
        G4double exp3p  = exp1+R22Apd;
        G4double NBinC  = N2p*BinCoeff;
        Din2           += NBinC*(R13Ap1/exp1p*std::exp(-Q24/exp1p)-
                                 R22Ap1/exp2p*std::exp(-Q24/exp2p)+
                                 R23Ap1/exp3p*std::exp(-Q24/exp3p));
        DmedTot        += NBinC * (R13Ap1/exp1p - R22Ap1/exp2p + R23Ap1/exp3p);
        N2p            *= -R23dR13;
      }
      G4double comFac   = N1p*Mnoj[i]/(i+2)/(i+1)*std::cos(FiH*i);
				  Din1             += Din2*comFac;
      DTot1            += DmedTot*comFac;
      if(std::abs(Din2*N1p/Din1) < eps) break;
    }
    G4double cFac       = A*(A-1)*4/Norm/Norm;
    Din1               *= -cFac;
    DTot1              *= cFac;
  }
  //else Din1=0.;

  G4double Corr0      = Tot00/Tot1;
  ImElA0             *= Corr0;
                        
  //return (ReElA0*ReElA0+(ImElA0+Din1)*(ImElA0+Din1))*CMMom*CMMom*mb2G2/4/pi2; // ds/do
  return (ReElA0*ReElA0+(ImElA0+Din1)*(ImElA0+Din1))*mb2G2/(twopi+twopi);
} // End of DiffElasticCS

G4double CHIPSDiffElCS(G4int hPDG, G4double HadrMom, G4int Z, G4int A, G4double aQ2) // MeV
{//      =================================================================
  static const G4double eps = 0.000001;                     // Accuracy of calculations
  static const G4double mb2G2 = 2.568;                      // Transform from mb to GeV^-2
  static const G4double piMG  = pi/mb2G2;                   // PiTrans from GeV^-2 to mb
  G4double p          = HadrMom/1000.;                      // Momentum (GeV, inputParam)
  G4double p2         = p*p;
		G4double MassH      = G4QPDGCode(hPDG).GetMass()/1000.;   // Hadron Mass in GeV
  G4double MassH2     = MassH*MassH;
  G4double HadrEnergy = std::sqrt(p2+MassH2);               // Tot energy in GeV
		//G4cout<<"G4GCS::DiffElasticCS: PGG_proj="<<hPDG<<", A_nuc= "<<A<<G4endl;
  //if(HadrEnergy-MassH<1.) G4Exception("Kin energy is too small for the model (T>1 GeV)");
  CalculateParameters(hPDG, HadrMom);                       // ?                   
  CalculateIntegralCS(A, HadrEnergy);                       // ?
  CalculateNuclearParameters(A);                            // ?
  G4double Q2 = aQ2/1000/1000;                              // in GeV^2
		//G4doubleMassN     = A*0.938;                              // @@ This is bad!
		G4double MassN    = A*0.9315;                             // @@ Atomic Unit is bad too
  if(G4NucleiPropertiesTable::IsInTable(Z,A))
				MassN=G4NucleiProperties::GetNuclearMass(A,Z)/1000.;    // Geant4 NuclearMass in GeV
  G4double MassN2   = MassN*MassN;
  G4double S        = (MassN+MassN)*HadrEnergy+MassN2+MassH2;// Mondelststam s
  G4double EcmH     = (S-MassN2+MassH2)/2/std::sqrt(S);     // CM energy of a hadron
  G4double CMMom    = std::sqrt(EcmH*EcmH-MassH2);          // CM momentum
  //G4cout<<"P="<<CMMom<<",E="<<HadrEnergy<<",N="<<MassN2<<",h="<<MassH2<<",p="<<p<<G4endl;
  G4double p3       = p2*p;
  G4double p4       = p2*p2;
  G4double sp       = std::sqrt(p);
  G4double p2s      = p2*sp;
  G4double ap       = std::log(p);
  G4double dl       = ap-3.;
  G4double dl2      = dl*dl;
  G4double shPPTot  = 2.91/(p2s+.0024)+5.+(32.+.3*dl2+23./p)/(1.+1.3/p4); // SigPP in mb
  G4double shPNTot  = 12./(p2s+.05*p+.0001/std::sqrt(sp))+.35/p
                      +(38.+.3*dl2+8./p)/(1.+1.2/p3);       // SigPP in mb
  G4double hNTot    = (Z*shPPTot+(A-Z)*shPNTot)/A;          // SighN in mb
  G4double Stot     = hNTot*mb2G2;                          // in GeV^-2
  G4double shPPSl   = 8.*std::pow(p,.055)/(1.+3.64/p4);     // PPslope in GeV^-2
  G4double shPNSl   = (7.2+4.32/(p4*p4+.012*p3))/(1.+2.5/p4); // PNslope in GeV^-2
  G4double Bhad     = (Z*shPPSl+(A-Z)*shPNSl)/A;            // B-slope in GeV^-2
  G4double hNReIm   = -.55+ap*(.12+ap*.0045);               // Re/Im_hN in no unit
  G4double Asq      = 1+hNReIm*hNReIm;                      // |M|^2/(ImM)^2
  G4double Rho2     = std::sqrt(Asq);                       // M/ImM
  G4double r1       = 3.9*std::pow(A-1.,.309);             // Positive diffractionalRadius
  G4double r2       = 2.*std::pow(A,.36);                  // Negative diffraction radius
  G4double pN       = Pnucl;                                // Dubna value
  //G4double pN       = .4;                                   // Screaning factor
  G4double Ae       = Aeff;                                 // Dubna value
  //G4double Ae       = .75;                                   // Normalization
  G4double R12      = r1*r1;
  G4double R22      = r2*r2;
  G4double R1C      = R12*r1;
  G4double R2C      = R22*r2;
  G4double R12B     = R12+Bhad+Bhad;                        // Slope is used
  G4double R22B     = R22+Bhad+Bhad;                        // Slope is used
  G4double Norm     = (R1C-pN*R2C)*Ae;                      // ScreanFac & NormFac are used
  G4double R13      = R1C/R12B;                             // Slope is used
  G4double R23      = pN*R2C/R22B;                          // Slope & ScreanFact are used
  G4double norFac   = Stot/twopi/Norm;                      // totCS (in GeV^-2) is used
  G4double Unucl    = norFac*R13;
  G4double SinFi    = hNReIm/Rho2;                          // Real part
  G4double FiH      = std::asin(SinFi);
  G4double N        = -1.;
  G4double N2       = R23/R13;                              // Slope is used
  G4double ImElA0   = 0.;
  G4double ReElA0   = 0.;
  G4double Tot1     = 0.;
  for(G4int i=1; i<=A; i++)
  {
    N              *= -Unucl*(A-i+1)*Rho2/i;                // TotCS is used
    G4double N4     = 1.;
    G4double Prod1  = std::exp(-Q2*R12B/i/4)*R12B/i;        // Slope is used
    G4double medTot = R12B/i;                               // Slope is used
    for(G4int l=1; l<=i; l++)
    {
      G4double exp1 = l/R22B+(i-l)/R12B;                    // Slope is used
      N4           *= -(i-l+1)*N2/l;                        // Slope is used
      G4double Nexp = N4/exp1;
      Prod1        += Nexp*std::exp(-Q2/exp1/4);
      medTot       += Nexp;
    }
    G4double FiHi   = FiH*i;
    G4double nCos   = N*std::cos(FiHi);
    ReElA0 += Prod1*N*std::sin(FiHi);
    ImElA0 += Prod1*nCos;
    Tot1           += medTot*nCos;
    if(std::abs(Prod1*N/ImElA0) < eps) break;
  }
  ImElA0         *= piMG;                          // The amplitude in mb
  ReElA0         *= piMG;                          // The amplitude in mb
  Tot1           *= piMG+piMG;
  G4double Corr0  = Tot00/Tot1;
  ImElA0         *= Corr0;
                        
  //return (ReElA0*ReElA0+ImElA0*ImElA0)*CMMom*CMMom*mb2G2/4/pi2; // ds/do
  return (ReElA0*ReElA0+ImElA0*ImElA0)*mb2G2/(twopi+twopi);
} // End of DiffElasticCS

G4double CHIPS_Tb(G4int  A, G4double b)              // T(b) in fm-2
{
  static G4int Am=0;                                 // Associative memory value A
  static G4double B=0.;                              // Effective edge
  static G4double C=0.;                              // Normalization factor
  static G4double D=0.;                              // Slope parameter
  if(A!=Am)
		{
    B=.0008*A*A;                                     // no units
    D=.42*std::pow(A,-.26);                               // fm^-2
    C=A*D/pi/std::log(1.+B);                         // fm^-2
  }
  G4double E=B*std::exp(-D*b*b);
  return C*E/(1.+E);                                 // T(b) in fm^-2
}

G4double Thickness(G4int A, G4double b, G4double R)  // T(b) in fm^-2
{//      ==========================================  // rAmax=2.5 is a GlobStatConstPar
  static const G4double bTh  = .545;                 // Surface thicknes
  static const G4double bTh2 = bTh*bTh;              // SquaredSurface thicknes
  G4double b2    = b*b;                              // Squared impact parameter
  G4double dr    = rAmax*R/(NptB-1);                 // Step of integration in z
  G4double R2    = R*R;                              // Working parameter
  G4double Norm  = .75/pi/R2/R/(1.+bTh2*pi2/R2);     // Corrected Volume (?)
  G4double r     = -dr;                              // Pre value for the radius
  G4double SumZ  = 0.;                               // Prototype of the result
  //G4double SumN  = 0.;                             // is not used
  for(G4int k=0; k<NptB; k++)
  {
    r    += dr;
    SumZ += 1./(1.+std::exp((std::sqrt(b2+r*r)-R)/bTh)); // Woods-Saxon density
    //SumN += r*r/(1.+std::exp((r-rAfm)/bTh));       // is not used
  }
  //SumN *= (twopi+twopi)*dr*Norm;                   // is not used
  return SumZ*(A+A)*dr*Norm;                         // T(b) in fm^-2
} // End of Thickness

G4double CoherentDifElasticCS(G4int hPDG, G4double HadrMom, G4int A, G4double aQ2) //in MeV
{//      =========================================================================
  //static const G4double eps = 0.000001;                   // CalculationAccuracy@@NotUsed
  //static const G4double mN = .9315;                      // Atomic Unit GeV
  static const G4double mN = .938;                      // Atomic Unit GeV
  static const G4double mb2G2 = 2.568;                      // Transform from mb to GeV^-2
  static const G4double incrf = 10;            // Increase factor of radius for intergation
  static const G4double m2G10 = incrf*mb2G2;                // Big conversion factor
  static const G4double s2G10 = std::sqrt(m2G10);           // sqrt Big conversion factor
  G4double ReIntegrand[NptB],ImIntegrand[NptB],Thick[NptB]; // Calculated arrays
  G4QBesIKJY QI0(BessI0);                                   // I0 Bessel function
  G4QBesIKJY QJ0(BessJ0);                                   // J0 Bessel function
  G4double hMom       = HadrMom/1000.;                      // Momentum (GeV, inputParam)
		G4double MassH      = G4QPDGCode(hPDG).GetMass()/1000.;   // Hadron Mass in GeV
  G4double MassH2     = MassH*MassH;
  G4double HadrEnergy = std::sqrt(hMom*hMom+MassH2);        // Tot energy in GeV
  if(A==2 || A==3) G4Exception("G4GCS: This model does not work for nuclei with A=2, A=3");
  if(A>252) G4Exception("G4GCS: This nucleus is too heavy for the model !!!");
  //if(HadrEnergy-MassH<1.) G4Exception("Kin energy is too small for the model (T>1 GeV)");
  G4double Re         = std::pow(A, 0.33333333);
  G4double Re2        = Re*Re;
  CalculateParameters(hPDG, HadrMom);                   
  G4double Q2 = aQ2/1000/1000;                              //  GeV
  // Individual corrections for nuclei which had data
  if    (A==208) r0   = 1.125;
  else if(A==90) r0   = 1.12;
  else if(A==64) r0   = 1.1;
  else if(A==58) r0   = 1.09;
  else if(A==48) r0   = 1.07;
  else if(A==40) r0   = 1.15;
  else if(A==28) r0   = 0.93;
  else if(A==16) r0   = 0.92;
  else if(A==12) r0   = 0.80;
		else           r0   = 1.16*(1.-1.16/Re2);  // For other nuclei which have not been tested
		G4double MassN      = A*0.9315;                             // @@ Atomic Unit is bad too
  G4double MassN2     = MassN*MassN;
  G4double S          = (MassN+MassN)*HadrEnergy+MassN2+MassH2; // Mondelstam S
  G4double EcmH       = (S-MassN2+MassH2)/2/std::sqrt(S);   // Hadron CM Energy
  G4double CMMom      = std::sqrt(EcmH*EcmH-MassH2);        // CM momentum
  G4double rAfm       = r0*Re;                              // Nuclear radius
  G4double rAGeV      = rAfm*s2G10;                         // Big integration radius
  G4double stepB      = rAmax*rAGeV/(NptB-1);               // dr step of integration
  G4double hTotG2     = HadrTot*mb2G2;
  G4double ReSum      = 0.;
  G4double ImSum      = 0.;
  G4double ValB       = -stepB;
  for(G4int i=0; i<NptB; i++)
  {
    ValB             += stepB;                              // An incident parameter
    G4double ValB2    = ValB*ValB;                          // A working value
    G4double IPH      = ValB/HadrSlope;                     // A working value
    G4double dHS      = HadrSlope+HadrSlope;                // A working value
    G4double Integ    = 0.;
    G4double ValS     = 0.;
    for(G4int j=1; j<NptB; j++)
    {
      ValS           += stepB;                              // Impact parameter GeV^-1
      if(!i) Thick[j] = Thickness(A,ValS/s2G10,rAfm)/m2G10; // Calculate only once
      G4double FunS   = IPH*ValS;
      if(FunS > 320.) break;
      Integ          += ValS*std::exp(-(ValS*ValS+ValB2)/dHS)*QI0(FunS)*Thick[j];
    } 
    G4double InExp    = -hTotG2*Integ*stepB/dHS;
    G4double expB     = std::exp(InExp);
    G4double HRIE     = HadrReIm*InExp;
    ReIntegrand[i]    = (1.-expB*std::cos(HRIE));           // Real part of the amplitude
    ImIntegrand[i]    = expB*std::sin(HRIE);                // Imaginary part of the amplit
  } 
  InCoh     = 0.;                                           // incohirent (quasi-elastic)
  ValB       = -stepB;
  for(G4int k=0; k<NptB; k++)                               // Third integration (?)
  {
    ValB              += stepB;
    InCoh             += Thick[k]*ValB*std::exp(-hTotG2*Thick[k]);
    G4double J0qb      = QJ0(std::sqrt(Q2)*ValB)*ValB;      // Bessel0(Q*b)
    ReSum             += J0qb*ReIntegrand[k];
    ImSum             += J0qb*ImIntegrand[k];
  }
  //G4cout<<"GHAD:st="<<stepB<<",hT="<<hTotG2<<",HRI="<<HadrReIm<<",HS="<<HadrSlope<<",Q2="
  //      <<Q2<<",m="<<mb2G2<<G4endl;
  InCoh *= stepB*hTotG2*hTotG2*(1.+HadrReIm*HadrReIm)*std::exp(-HadrSlope*Q2)/8/mb2G2;
  //G4cout<<"GHAD:RS="<<ReSum<<",IS="<<ImSum<<",CM="<<CMMom<<",st="<<stepB<<G4endl;
  //return (ReSum*ReSum+ImSum*ImSum)*m2G10*CMMom*CMMom*stepB*stepB/twopi; // ds/do
  return (ReSum*ReSum+ImSum*ImSum)*m2G10*stepB*stepB/12; // ds/dt
}  

G4double CHIPSDifElasticCS(G4int hPDG, G4double hMom, G4int A, G4int Z, G4double aQ2)
{//      ============================================================================
  static const G4int Npb      = 500;                        // A#of intergation points
  //static const G4double mN = .9315;                         // Atomic Unit GeV
  static const G4double mN = .938;                          // Atomic Unit GeV
  static const G4double hc2   = .3893793;                   // Transform from GeV^-2 to mb
  static const G4double mb2G2 = 1./hc2;                     // Transform from mb to GeV^-2
  static const G4double f22mb = 10;                         // Transform from fermi^2 to mb
  static const G4double f22G2 = f22mb*mb2G2;                // Transform from fm2 to GeV^-2
  static const G4double f2Gm1 = std::sqrt(f22G2);           // Transform from fm to GeV^-1
  G4double Re         = std::pow(A,.33333333);              // A-dep coefficient
  G4double Lim        = 100*Re;                             // Integration accuracy limit
  G4double Tb[Npb];                                         // Calculated T(b) array
  G4QBesIKJY QI0(BessI0);                                   // I0 Bessel function
  G4QBesIKJY QJ0(BessJ0);                                   // J0 Bessel function
  G4double p          = hMom/1000.;                         // Momentum (GeV, inputParam)
  G4double p2         = p*p;
		G4double MassH      = G4QPDGCode(hPDG).GetMass()/1000.;   // Hadron Mass in GeV
  G4double MassH2     = MassH*MassH;                        // Squared mass of the hadron
  G4double hEnergy    = std::sqrt(p2+MassH2);               // Tot energy in GeV
  G4double Q2         = aQ2/1000000.;                       // -t in GeV
		G4double MassN    = A*0.9315;                             // @@ Atomic Unit is bad too
  if(G4NucleiPropertiesTable::IsInTable(Z,A))
				MassN=G4NucleiProperties::GetNuclearMass(A,Z)/1000.;    // Geant4 NuclearMass in GeV
  G4double MassN2     = MassN*MassN;                        // Squared mass of the target
  G4double S          = (MassN+MassN)*hEnergy+MassN2+MassH2;// Mondelstam s
  G4double EcmH       = (S-MassN2+MassH2)/2/std::sqrt(S);   // Hadron CM Energy
  G4double CMMom      = std::sqrt(EcmH*EcmH-MassH2);        // CM momentum (to norm CS)
  //G4cout<<"2: P="<<CMMom<<",E="<<hEnergy<<",N="<<MassN2<<",h="<<MassH2<<",p="<<p<<G4endl;
  // The mean value of the total can be used
  G4double p3         = p2*p;
  G4double p4         = p2*p2;
  G4double sp         = std::sqrt(p);
  G4double p2s        = p2*sp;
  G4double ap         = std::log(p);
  G4double dl         = ap-3.;
  G4double dl2        = dl*dl;
  G4double shPPTot    = 2.91/(p2s+.0024)+5.+(32.+.3*dl2+23./p)/(1.+1.3/p4); // SigPP in mb
  G4double shPNTot    = 12./(p2s+.05*p+.0001/std::sqrt(sp))+.35/p
                        +(38.+.3*dl2+8./p)/(1.+1.2/p3);     // SigPP in mb
  G4double shNTot     = (Z*shPPTot+(A-Z)*shPNTot)/A;        // SighN in mb
#ifdef debug
  G4cout<<"CHIPS:SI,p="<<p<<",n="<<shNTot<<",P="<<shPPTot<<",N="<<shPNTot
        <<",Z="<<Z<<",A="<<A<<G4endl;
#endif
  G4double shPPSl     = 8.*std::pow(p,.055)/(1.+3.64/p4);   // PPslope in GeV^-2
  G4double shPNSl     = (7.2+4.32/(p4*p4+.012*p3))/(1.+2.5/p4); // PNslope in GeV^-2
  G4double shNSl      = (Z*shPPSl+(A-Z)*shPNSl)/A;          // B-slope in GeV^-2
#ifdef debug
  G4cout<<"CHIPS:SL,n="<<shNSl<<",P="<<shPPSl<<",N="<<shPNSl<<G4endl;
#endif
  G4double shNReIm    = -.55+ap*(.12+ap*.0045); // Re/Im_hN in no unit
  //G4cout<<"CHPS: s="<<S<<",T="<<shNTot<<",R="<<shNReIm<<",B="<<shNSl<<G4endl;
  // @@ End of temporary ^^^^^^^
  G4double dHS        = shNSl+shNSl;                        // Working: doubled B-slope
  G4double stepB      = (Re+Re+2.7)*f2Gm1/(Npb-1);          // in GeV^-1, step of integral
  G4double hTotG2     = shNTot*mb2G2;                       // sigma_hN in GeV^-2
  G4double ReSum      = 0.;                                 // Integration of RePart of Amp
  G4double ImSum      = 0.;                                 // Integration of ImPart of Amp
  G4double ValB       = -stepB;
  for(G4int i=0; i<Npb; i++)                                // First integration over b
  {
    ValB             += stepB;                              // An incident parameter
    G4double ValB2    = ValB*ValB;                          // A working value b^2
    G4double IPH      = ValB/shNSl;                         // A working value slope
    G4double Integ    = 0.;                                 // Integral over ImpactParam.
    G4double ValS     = 0.;                                 // Prototype of ImpactParameter
    for(G4int j=1; j<Npb; j++)                              // Second integration over b
    {
      ValS           += stepB; //  back to fm               // Impact parameter GeV^-1
      if(!i) Tb[j]    = CHIPS_Tb(A,ValS/f2Gm1)/f22G2;       // GeV^2, calculate only once
      //if(!i) Tb[j]    = Thickness(A,ValS/f2Gm1,rAfm)/f22G2; // Calculate T(b) only once
      G4double FunS   = IPH*ValS;                           // b1*b2/slope
      if(FunS > Lim) break;                                 // To avoid NAN
      Integ          += ValS*std::exp(-(ValS*ValS+ValB2)/dHS)*QI0(FunS)*Tb[j]; // BessI0
    } 
    G4double InExp    = -hTotG2*Integ*stepB/dHS;            // Integrated absorption 
    G4double expB     = std::exp(InExp);                    // Exponential absorption
    G4double HRIE     = shNReIm*InExp;                      // Phase shift
    G4double J0qb      = QJ0(std::sqrt(Q2)*ValB)*ValB;      // Bessel0(Q*b)
    ReSum             += J0qb*(1.-expB*std::cos(HRIE));
    ImSum             += J0qb*expB*std::sin(HRIE);
  } 
  //G4cout<<"CHPS:RS="<<ReSum<<",IS="<<ImSum<<",CM="<<CMMom<<",st="<<stepB<<G4endl;
  //return (ReSum*ReSum+ImSum*ImSum)*mb2G2*CMMom*CMMom*stepB*stepB/twopi;
  //return (ReSum*ReSum+ImSum*ImSum)*f22G2*CMMom*CMMom*stepB*stepB/twopi; // ds/do
  return (ReSum*ReSum+ImSum*ImSum)*f22G2*stepB*stepB/12; // ds/dt
}

// Separate quasielastic calculation
G4double CHIPSDifQuasiElasticCS(G4int hPDG, G4double hMom, G4int A, G4int Z, G4double aQ2)
{//      =================================================================================
  static const G4int Npb      = 500;                        // A#of intergation points
  static const G4double hc2   = .3893793;                   // Transform from GeV^-2 to mb
  static const G4double mb2G2 = 1./hc2;                     // Transform from mb to GeV^-2
  //static const G4double mb2G2 = 2.568;                     // Transform from mb to GeV^-2
  static const G4double f22mb = 10;                         // Transform from fermi^2 to mb
  static const G4double f22G2 = f22mb*mb2G2;                // Transform from fm2 to GeV^-2
  static const G4double f2Gm1 = std::sqrt(f22G2);           // Transform from fm to GeV^-1
  G4double p          = hMom/1000.;                         // Momentum (GeV, inputParam)
  G4double p2         = p*p;
  G4double Q2         = aQ2/1000000.;                       // -t in GeV
  // @@ Temporary only for nucleons
  G4double p3         = p2*p;
  G4double p4         = p2*p2;
  G4double sp         = std::sqrt(p);
  G4double p2s        = p2*sp;
  G4double ap         = std::log(p);
  G4double dl         = ap-3.;
  G4double dl2        = dl*dl;
  G4double shPPTot    = 2.91/(p2s+.0024)+5.+(32.+.3*dl2+23./p)/(1.+1.3/p4); // SigPP in mb
  G4double shPNTot    = 12./(p2s+.05*p+.0001/std::sqrt(sp))+.35/p
                        +(38.+.3*dl2+8./p)/(1.+1.2/p3);     // SigPP in mb
  G4double shNTot     = (Z*shPPTot+(A-Z)*shPNTot)/A;        // SighN in mb
  G4double shPPSl     = 8.*std::pow(p,.055)/(1.+3.64/p4);   // PPslope in GeV^-2
  G4double shPNSl     = (7.2+4.32/(p4*p4+.012*p3))/(1.+2.5/p4); // PNslope in GeV^-2
  G4double shNSl      = (Z*shPPSl+(A-Z)*shPNSl)/A;          // B-slope in GeV^-2
  G4double shNReIm    = -.55+ap*(.12+ap*.0045); // Re/Im_hN in no unit
  // @@ End of temporary ^^^^^^^
  G4double Re         = std::pow(A,.33333333);              // A-dep coefficient
  G4double stepB      = (Re+Re+2.7)*f2Gm1/(Npb-1);          // in GeV^-1, step of integral
  G4double InQE       = 0.;                                 // Quasielastic integral
  G4double ValB       = -stepB;                             // Impact parameter
  G4double hTotG2     = shNTot*mb2G2;                       // sigma_hN in GeV^-2
  for(G4int k=0; k<Npb; k++)
  {
    ValB        += stepB;
    //G4double Tb  = Thickness(A,ValB/f2Gm1,rAfm)/f22G2;      // Calculate T(b) only once
    G4double Tb  = CHIPS_Tb(A,ValB/f2Gm1)/f22G2;            // GeV^2, calculate only once
    InQE        += Tb*ValB*std::exp(-hTotG2*Tb);
  }
  //G4cout<<"CHPS:st="<<stepB<<",hT="<<hTotG2<<",HRI="<<HadrReIm<<",HS="<<HadrSlope<<",Q2="
  //      <<Q2<<",m="<<mb2G2<<G4endl;
  return InQE*stepB*hTotG2*hTotG2*(1.+shNReIm*shNReIm)*std::exp(-shNSl*Q2)/8/mb2G2;
}
#endif

int main()
{
#ifdef bestest
  // *** Test of CHIPS implementation of Bessel functions (accuracy 1.E-8 and better) ***
  const G4int nj=21;
  const G4int ny=16;
  const G4int ni=20;
  const G4int nk=20;
  const G4double jX[nj]=
          				{-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};
  const G4double yX[ny]={.1,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};
  const G4double iX[ni]=
                    {0.,.2,.4,.6,.8,1.,1.2,1.4,1.6,1.8,2.,2.5,3.,3.5,4.,4.5,5.,6.,8.,10.};
  const G4double kX[nk]=
                    {.1,.2,.4,.6,.8,1.,1.2,1.4,1.6,1.8,2.,2.5,3.,3.5,4.,4.5,5.,6.,8.,10.};
  const G4double vBJ0[nj]=
				              {-.1775968,-.3971498,-.2600520,0.2238908,0.7651976,1.0000000,0.7651977,
                   0.2238908,-.2600520,-.3971498,-.1775968,0.1506453,0.3000793,0.1716508,
                   -.0903336,-.2459358,-.1711903,0.0476893,0.2069261,0.1710735,-.0142245};
  const G4double vBJ1[nj]=
              				{0.3275791,0.0660433,-.3390590,-.5767248,-.4400506,0.0000000,0.4400506,
                   0.5767248,0.3390590,-.0660433,-.3275791,-.2766839,-.0046828,0.2346364,
                   0.2453118,0.0434728,-.1767853,-.2234471,-.0703181,0.1333752,0.2051040};
  const G4double vBY0[ny]=
        {-1.534239,0.0882570,.51037567,.37685001,-.0169407,-.3085176,-.2881947,-.0259497,
         0.2235215,0.2499367,0.0556712,-.1688473,-.2252373,-.0782079,0.1271926,0.2054743};
  const G4double vBY1[ny]=
    				{-6.458951,-.7812128,-.1070324,0.3246744,0.3979257,0.1478631,-.1750103,-.3026672,
         -.1580605,0.1043146,0.2490154,0.1637055,-.0570992,-.2100814,-.1666448,0.0210736};
  const G4double vBI0[ni]={1.0000000,1.0100250,1.0404018,1.0920453,1.1665149,
                           1.2660658,1.3937256,1.5533951,1.7499807,1.9895593,
                           2.2795852,3.2898391,4.8807925,7.3782035,11.301922,
                           17.481172,27.239871,67.234406,427.56411,2815.7167};
  const G4double vBI1[ni]={.00000000,.10050083,.20402675,.31370403,.43286480,
                           .56515912,.71467794,.88609197,1.0848107,1.3171674,
                           1.5906369,2.5167163,3.9533700,6.2058350,9.7594652,
                           15.389221,24.335643,61.341937,399.87313,2670.9883};
  const G4double vBK0[nk]={2.4270690,1.7527038,1.1145291,.77752208,.56534710,
                           .42102445,.31850821,.24365506,.18795475,.14593140,
                           .11389387,.062347553,.0347395,.019598897,.011159676,
                           .0063998572,.0036910983,.0012439943,.00014647071,.000017780062};
  const G4double vBK1[nk]={9.8538451,4.7759725,2.1843544,1.3028349,.86178163,
                           .60190724,.43459241,.32083589,.24063392,.18262309,
																											.13986588,.073890816,.040156431,.022239393,.012483499,
                           .0070780949,.0040446134,.0013439197,.00015536921,.000018648773};
  //Bessel J0('J',0); // CLHEP J/Y is possible insead of 'J'
  //Bessel J1('J',1);
  G4QBesIKJY QI0(BessI0); // CHIPS I/K/J/Y 0/1 Bessel functions
  G4QBesIKJY QI1(BessI1);
  G4QBesIKJY QJ0(BessJ0);
  G4QBesIKJY QJ1(BessJ1);
  G4QBesIKJY QK0(BessK0);
  G4QBesIKJY QK1(BessK1);
  G4QBesIKJY QY0(BessY0);
  G4QBesIKJY QY1(BessY1);
  G4cout<<"G4GCS: ***J0*** Test: ================================================"<<G4endl;
  for(G4int j0=0; j0<nj; j0++)
		{
    G4double F=QJ0(jX[j0]);
    G4double d=F-vBJ0[j0];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<jX[j0]<<", J0="<<F<<" - "<<vBJ0[j0]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***J1*** Test: ================================================"<<G4endl;
  for(G4int j1=0; j1<nj; j1++)
		{
    G4double F=QJ1(jX[j1]);
    G4double d=F-vBJ1[j1];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<jX[j1]<<", J1="<<F<<" - "<<vBJ1[j1]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***Y0*** Test: ================================================"<<G4endl;
  for(G4int y0=0; y0<ny; y0++)
		{
    G4double F=QY0(yX[y0]);
    G4double d=F-vBY0[y0];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<yX[y0]<<", Y0="<<F<<" - "<<vBY0[y0]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***Y1*** Test: ================================================"<<G4endl;
  for(G4int y1=0; y1<ny; y1++)
		{
    G4double F=QY1(yX[y1]);
    G4double d=F-vBY1[y1];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<yX[y1]<<", Y1="<<F<<" - "<<vBY1[y1]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***I0*** Test: ================================================"<<G4endl;
  for(G4int i0=0; i0<ni; i0++)
		{
    G4double F=QI0(iX[i0]);
    G4double d=F-vBI0[i0];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<iX[i0]<<", I0="<<F<<" - "<<vBI0[i0]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***I1*** Test: ================================================"<<G4endl;
  for(G4int i1=0; i1<ni; i1++)
		{
    G4double F=QI1(iX[i1]);
    G4double d=F-vBI1[i1];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<iX[i1]<<", I1="<<F<<" - "<<vBI1[i1]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***K0*** Test: ================================================"<<G4endl;
  for(G4int k0=0; k0<nk; k0++)
		{
    G4double F=QK0(kX[k0]);
    G4double d=F-vBK0[k0];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<kX[k0]<<", I1="<<F<<" - "<<vBK0[k0]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"G4GCS: ***K1*** Test: ================================================"<<G4endl;
  for(G4int k1=0; k1<nk; k1++)
		{
    G4double F=QK1(kX[k1]);
    G4double d=F-vBK1[k1];
    G4double r=std::abs(d/F);
    G4cout<<"G4GCS: x="<<kX[k1]<<", I1="<<F<<" - "<<vBK1[k1]<<" = "<<d<<", r="<<r<<G4endl;
  }
  G4cout<<"End .................................................................."<<G4endl;
#endif

  // Test of integrated cross sections
  //const G4int na=12;
  const G4int na=1;
  //
  const G4int np=7;
  const G4int nm=1;
  ////                He Be C O  Al Ti Ni Cu  Sn Ta  Pb   U
  //const G4int A[na]={4,9,12,16,27,48,58,64,120,181,207,238}; // A's of target nuclei
  //                 He Al Pb
  const G4int A[na]={119}; // A's of target nuclei
  //                     p    n  pi+  pi-  K+  K-  antip
  const G4int pdg[np]={2212,2112,211,-211,321,-321,-2212}; // projectiles
  const G4double mom[nm]={120.}; // momentum in MeV/c
#ifdef integrc
  for(G4int ip=0; ip<np; ip++)
		{
    G4int PDG=pdg[ip];
    for(G4int im=0; im<nm; im++)
		  {
      CalculateParameters(PDG, mom[im]);
      G4cout<<G4endl<<"G4GCS:PDG="<<PDG<<",P="<<mom[im]<<":A,Tot,Inel,Prod,El,QEl"<<G4endl;
      for(G4int ia=0; ia<na; ia++)
		    {
        CalculateIntegralCrossSections(A[ia]);
        G4cout<<A[ia]<<" "<<TotalCrSec<<" "<<InelCrSec<<" "<<ProdCrSec
              <<" "<<ElasticCrSec<<" "<<QuasyElasticCrSec<<G4endl;
      }
    }
  }
#endif

#ifdef eldiffcs
  const G4double ms=.000001;
  const G4int nt=200;
  G4double t[nt]; // Q2=-t in Mev^2
  const G4double tMin=200.;
  const G4double tMax=2000000.;
  const G4double ltMi=std::log(tMin);
  const G4double ltMa=std::log(tMax);
  const G4double dlt=(ltMa-ltMi)/(nt-1);
  G4double lt=ltMi-dlt;
  for(G4int it=0; it<nt; it++)
		{
    lt+=dlt;
    t[it]=std::exp(lt);
  }
  ////                He Be C O  Al Ti Ni Cu  Sn Ta  Pb   U
  //const G4int Z[na]={2,4, 6, 8,13,22,28,29, 50, 73, 82, 92}; // Z's of target nuclei
  //                He Al Pb
  const G4int Z[na]={50}; // Z's of target nuclei
  // Test of differential ellastic cross sections
		Initialize();
  //for(G4int ip=0; ip<np; ip++)
  for(G4int ip=0; ip<1; ip++)
		{
    G4int PDG=pdg[ip];
    for(G4int im=0; im<nm; im++)
		  {
      CalculateParameters(PDG, mom[im]);
      G4cout<<G4endl<<"G4GCS:PDG="<<PDG<<",P="<<mom[im]<<G4endl;
      for(G4int ia=0; ia<na; ia++)
		    {
        G4cout<<"G4GCS:-------------------A="<<A[ia]<<G4endl;
        for(G4int it=0; it<nt; it++)
		      {
          G4double Sig1 = CHIPSDiffElCS(PDG, mom[im], Z[ia], A[ia], t[it]);
          //G4double Sig1 = DiffElasticCS(PDG, mom[im], A[ia], t[it]); //Doesn't work
          G4double Sig2 = CoherentDifElasticCS(PDG, mom[im], A[ia], t[it]);
          G4double Sig3 = CHIPSDifElasticCS(PDG, mom[im], A[ia], Z[ia], t[it]);
          G4double CQEl = CHIPSDifQuasiElasticCS(PDG, mom[im], A[ia], Z[ia], t[it]);
								  G4cout<<std::setw(12)<<ms*t[it]<<std::setw(12)<<Sig1<<std::setw(12)<<Sig2
                <<std::setw(12)<<InCoh<<std::setw(12)<<Sig3<<std::setw(12)<<CQEl<<G4endl;
        }
      }
    }
  }
#endif
  return 0;
}
