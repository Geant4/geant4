//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// G4 Tools program: PhotoNuclearCalculation of gamma+A cross sections
// according to the paper M.V.Kossov "Approximation of photonuclear
// interaction cross sections" (EPJA,2002, be published). The cross
// sections are calculated according to the exact complicated formulas
// and tabulated for the fast response. The total energy range is
// subdivided in three regions:
// a) The GDR region (from the hadron production threshold to 106 MeV)
//    covering 46 nuclei: H2,He4,Li6,Li7,Be9,C12,N14,N15,O16,F19,Na23,
//    Mg24,Al27,Si28,S32,S34,(Ar40,Ca40),Fe54,Mn55,Fe56,Ni58,Co59,Cu,Zn
//    Se76,Se82,Ag,Cd,Sn,I,Sm154,Gd156,Tb159,Ho165,Er168,Yb174,Hf178,Hf180,
//    Ta181,W184,W186,Au197,Tl,Pb,Bi209,Th232,U235,U238,Pu239. For these
//    isotops the approximation fits the existing measurements. For
//    them the basic functions are calculated. For other isotops the
//    approximation is interpolated linearly with A.
// b) The resonance region (from 106 MeV to 50 GeV)
//    covering 14 nuclei: H1,H2,He3,He4,Li6,Li7,Be,C,O,Al,Cu,Sn,Pb,U.
//    For these isotops the melting of resonances in nuclear matter is
//    measured and approximated. For them the basic functions are
//    calculated. For other isotops the approximation is interpolated
//    linearly with A.
// c) The high energy region is an A-dependent function for E>100 GeV,
//    taking into account nuclear shaddowing.
// .....................................................
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-May-02
// The lust update: M.V. Kossov, CERN/ITEP(Moscow) 17-June-02
// 
//=====================================================================
#include "globals.hh"
#include "g4std/iostream"
#include "g4std/fstream"
#include "g4std/vector"
#include "G4ios.hh"

int main()
{
  const G4int mN=49;                // A#of GDR basic nuclei
  const G4int mC=105;               // A#of Resonance points in E (each MeV from 2 to 106)
  const G4double mA[mN]={
	2.,4.,6.,7.,9.,12.,14.,15.,16.,19.,23.,24.,27.,28.,32.,34.,40.,54.,
    55.,56.,58.7,58.9,63.5,65.4,76.,82.,107.9,112.4,118.7,126.9,154.,156.,159.,165.,
    168.,174.,178.,180.,181.,184.,186.,197.,204.4,207.2,209.,232.,235.,
    238.,239.};
  const G4int nN=14;                // A#of Resonance basic nuclei
  const G4int nC=224;               // A#of Resonance points in lnE
  const G4double nA[nN]={1.,2.,3.,4.,6.,7.,9.,12.,16.,27.,63.5,118.7,207.2,238.};
  const G4double r4[6]={0.,3.38,3.51,3.23,3.57,3.60};
  const G4double t4[6]={0.,3.14,3.09,2.80,3.25,3.17};
  const G4double r8[mN]={
    0.00,0.00,0.00,2.57,3.26,3.52,3.51,3.51,3.49,3.05,3.20,3.25,3.40,3.42,3.38,3.50,3.51,
    3.59,3.45,3.46,3.52,3.56,3.40,3.40,3.29,3.38,3.50,3.49,3.41,3.46,3.40,3.37,3.40,3.29,
    3.37,3.21,3.45,3.45,3.39,3.36,3.36,3.39,3.47,3.40,3.36,3.31,3.30,3.37,3.36};
  const G4double t8[mN]={
    0.00,0.00,0.00,2.48,3.05,3.11,3.09,3.08,3.10,2.70,2.90,2.95,3.00,2.97,2.98,3.00,3.05,
    2.94,2.93,2.86,2.99,3.01,2.85,2.84,2.72,2.76,2.81,2.77,2.72,2.75,2.78,2.74,2.72,2.72,
    2.70,2.74,2.70,2.70,2.64,2.63,2.61,2.59,2.63,2.57,2.63,2.62,2.62,2.67,2.65};
  const G4double iE=log(106.);      // Start logarithm energy
  const G4double fE=log(50000.);    // Finish logarithm energy (each 2.75 percent)
  const G4double dE=(fE-iE)/(nC-1); // Step in logarithm energy
  G4std::ofstream fileGDR("GDR.out", G4std::ios::out);
  fileGDR.setf( G4std::ios::scientific, G4std::ios::floatfield );
  G4std::ofstream fileRes("Res.out", G4std::ios::out);
  fileRes.setf( G4std::ios::scientific, G4std::ios::floatfield );
  G4int np=0;
  for(G4int m=0; m<mN; m++)
  {
    G4double A=mA[m];
    G4double lnA=log(A);
    G4double A2=A*A;
    G4double red=1.+16/A2/A2;
    G4double rho1=(3.2+.75*lnA)/red;
    G4double tau1=(6.6-.5*lnA)/red;
	if(A==2.)
	{
	  rho1=1.86;
	  tau1=1.2;
	}
	G4double rho2=(4.+.125*lnA)/red;
	G4double tau2=3.4/red;
	if(A==6.)
	{
	  rho2=2.9;
	  tau2=2.32;
	}
	else if(A==2.)
	{
	  rho2=2.11;
	  tau2=1.5;
	}
	G4double rho4=3.8+.05*lnA;
	G4double tau4=3.8-.25*lnA;
	if(A<13.)
	{
	  rho4=r4[m];
	  tau4=t4[m];
	}
	G4double rho8=r8[m];
	G4double tau8=t8[m];
    G4double tr=5.13-.00075*A; // Resonance threshold Position
    //if(A==1.)tr=5.24;
	G4double rr=11.;           // Resonance threshold Slope
	if(A<2.5)rr=25.;
	G4double dw=.056+lnA*(.03-.001*lnA); // Delta width
	G4double du=5.82-.07/(1.+.003*A2);   // Delta position
    G4double da=.39*A;                   // Delta amplitude
    if(A==2.)da=.88;
    //else if(A==1.) da=.55;
    G4double ekin=1.;
    fileGDR<<"  static const G4double SL"<<m<<"[nL]={"<<G4endl<<"    ";
    G4cout<<"**** A_low="<<A<<G4endl;
    np=0;
	for(G4int em=0; em<mC; em++)
	{
      ekin+=1.;
      G4double z=log(ekin);
      G4double ds=z-du;
	  G4double sigd=da/(1.+ds*ds/dw)/(1.+exp(rr*(tr-z))); // Delta contribution
	  G4double g1=exp(rho1-z)/(1.+exp(3*(tau1-z)));
	  G4double g2=exp(2*(rho2-z))/(1.+exp(6*(tau2-z)));
      G4double g4=0.;
      if(A>3.5)g4=exp(4*(rho4-z))/(1.+exp(12*(tau4-z)));
      G4double g8=0.;
      if(A>6.5)g8=exp(8*(rho8-z))/(1.+exp(24*(tau8-z)));
      G4double sig=sigd+g1+g2+g4+g8;
      np++;
      if(np==7)
      {
        if(em==mC-1) fileGDR<<sig<<"};"<<G4endl;
        else         fileGDR<<sig<<","<<G4endl<<"    ";
	  }
      else      fileGDR<<sig<<",";
      if(np==7) np=0;
    } // End of the point LOOP
  } // End of the isotop LOOP
  // ===== High Energy Calculations =======================
  //G4double shd=1.0663-.0023*log(2.);   // HE PomShadowing(D)
  np=0;
  for(G4int n=0; n<nN; n++)
  {
    G4double A=nA[n];
    G4double A2=A*A;
    fileRes<<"  static const G4double SH"<<n<<"[nH]={"<<G4endl<<"    ";
    G4cout<<"**** A_high="<<A<<G4endl;
    G4double lnA=log(A);
    G4double slA=sqrt(lnA);
    G4double red=1.+16/A2/A2;
    G4double rho1=(3.2+.75*lnA)/red;
    G4double tau1=(6.6-.5*lnA)/red;
	if(A==2.)
	{
	  rho1=1.86;
	  tau1=1.2;
	}
	G4double rho2=(4.+.125*lnA)/red;
	G4double tau2=3.4/red;
	if(A==6.)
	{
	  rho2=2.9;
	  tau2=2.32;
	}
	else if(A==2.)
	{
	  rho2=2.11;
	  tau2=1.5;
	}
    G4double rho4=6.27;
    G4double tau4=7.25;
    G4double rho8=6.66;
    G4double tau8=6.90;
    if(A==2.)
	{
     rho4=6.2;
     tau4=7.1;
     rho8=6.62;
     tau8=6.91;
    }
    G4double tr=5.13-.00075*A; // Resonance threshold Position
    if(A==1.)tr=5.24;
	G4double rr=11.;           // Resonance threshold Slope
	if(A<2.5)rr=25.;
	// Delta
	G4double dw=.056+lnA*(.03-.001*lnA); // Delta width
	G4double du=5.82-.07/(1.+.003*A2);   // Delta position
    G4double da=.39*A;                   // Delta amplitude
    if(A==2.)da=.88;
    else if(A==1.) da=.55;
	// High Resonance
	G4double hw=.045+.04*slA*slA*slA;    // HighR width
	G4double hu=6.496+.042*lnA;          // HighR position
    if(A==1.)hu=6.57;
    else if(A==2.)hu=6.575;
    G4double ha=.223;                    // HighR amplitude
    if(A>2.5)ha=.16*A/slA;
    else if(A==2.) ha=.348;
    G4double sp=A*(1.-.072*lnA);         // HE TotShadowing
    G4double sh=1.0663-.0023*lnA;        // HE PomShadowing
    if(A==1.)sh=1.07;
	// -- Loop over energies ---
    G4double z=iE-dE;
    np=0;
	for(G4int en=0; en<nC; en++)
	{
      z+=dE;
      G4double ds=z-du;
      G4double fr=1.+exp(rr*(tr-z));
	  G4double sigd=da/(1.+ds*ds/dw); // Delta contribution
      G4double hs=z-hu;
	  G4double sigh=ha/(1.+hs*hs/hw); // HighR contribution
	  G4double g1=0.;
      if(A>1.5)g1=exp(rho1-z)/(1.+exp(3*(tau1-z)));
	  G4double g2=0.;
      if(A>1.5)g2=exp(2*(rho2-z))/(1.+exp(6*(tau2-z)));
      G4double g4=0.;
      if(A<2.5)g4=exp(4*(rho4-z))/(1.+exp(12*(tau4-z)));
      G4double g8=0.;
      if(A<2.5)g8=exp(8*(rho8-z))/(1.+exp(24*(tau8-z)));
      G4double hp=.0375*(z-16.5)+sh*exp(-.11*z);
      G4double fp=hp/(1.+exp(4*(7.-z)));
      G4double sig=(sigd+sigh)/fr+g1+g2+g4+g8+sp*fp;
      np++;
      if(np==7)
      {
        if(en==nC-1) fileRes<<sig<<"};"<<G4endl;
        else         fileRes<<sig<<","<<G4endl<<"    ";
	  }
      else      fileRes<<sig<<",";
      if(np==7) np=0;
      //if(en==nC)
	  //{
      //  G4double fsig=sp*(.0375*(z-16.5)+shd*exp(-.11*z));
      //  fileRes<<">>> fun="<<fsig<<G4endl;
	  //}
    } // End of the point LOOP
  } // End of the isotop LOOP
  return EXIT_SUCCESS;
}
