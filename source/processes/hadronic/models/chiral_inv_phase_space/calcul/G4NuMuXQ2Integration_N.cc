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
// G4 Tools program: NuMu DIS (x,Q2) approximation is integrated over x & Q2
// .....................................................[ for a <nucleon> ]
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Oct-07
// 
//=====================================================================
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include "G4ios.hh"
#include <CLHEP/GenericFunctions/LogGamma.hh>

void strucf(int A, double x, double Q2, double& f2, double& xf3, double& fL)
{
  //const double MN=.931494043;   // Nucleon mass (inside nucleus, atomic mass unit, GeV)
  //const double MN2=MN*MN;       // M_N^2 in GeV^2
  //const double mpi=.13957018;   // charged pi meson mass in GeV
  //const double Wt=MN+mpi;       // Delta threshold
  //const double W2t=Wt*Wt;       // Squared Delta threshold
  const Genfun::LogGamma lGam;
  static int mA=0;
  static double mQ2=0., mN, mD, mDel, mU2, mU3, mV;
  //static double mUU;
  double N=3., D=0., Del=0., U2=0., U3=0., V=0.;
  //double UU=0.;
  if(A==mA && Q2==mQ2)         // Associative memory for acceleration
		{
    N  =mN;
    D  =mD;
    Del=mDel;
    U2 =mU2;
    U3 =mU3;
    //UU =mUU;
    V  =mV;
  }
  else
  {
    double Q=std::sqrt(Q2);
    N=3.+.3581*std::log(1.+Q2/.04); // a#of partons in the nonperturbative phase space
    Del=.077+.393/(1.+5./Q);
    D=.68*std::pow(1.+.145/Q2,-1.-Del);
    V=3*(N-1.);
    double c3=.75;
    double uu=std::exp(lGam(N-Del)-lGam(N-1.)-lGam(1.-Del))/N;
    U2=(c3+N-3.)*uu;
    U3=c3*uu;
    //UU=uu+uu+uu; // @@
    mA  = A;
    mQ2 = Q2;
    mN  = N;
    mD  = D;
    mDel=Del;
    mU2 =U2;
    mU3 =U3;
    //mUU =UU; // @@
    mV  =V;
  }  
  // From here the Q2 coefficients are used
  double x1=std::pow(1.-x,N-2.);
  double pp=D*std::pow(x,-Del)*x1;
  double dir=(1-D)*V*x*x1;
  double per=U2*pp;
		f2 = per + dir;
  //double W2=MN2-MN2*x+Q2/x-Q2;
  //if(W2<W2t)
		//{
  //  per=UU*pp;
  //  xf3= per+dir;
  //}
  //else
  xf3= U3*pp+dir;
  fL = per/5.;
  return;
}

void getFun(int A, double lx, double Q2, double* f)
{
  const double MN=.931494043;   // Nucleon mass (inside nucleus, atomic mass unit, GeV)
  const double MN2=MN*MN;   // Squared Nucleon mass (inside nucleus, atomic mass unit, GeV)
  double f2=0., xf3=0., fL=0.;
  if (lx>0.5) G4cerr<<"***getFun: ln(x)="<<lx<<">.5"<<G4endl;
  double x=std::exp(lx);
  double x2=x*x;
  double c=0.;
  if(Q2>.0001) c=x2*MN2/Q2;
  double c2=1+c+c;     // Q^2/nu^2 correction
  strucf(A, x, Q2, f2, xf3, fL);
  f[0]=f2;             // direct part
  f[1]=(-f2+xf3)/x;    // *y (neutrino) part
  f[2]=(-f2-xf3)/x;    // *y (anti-neutrino) part
  f[3]=(f2*c2-fL-xf3)/x2; // *y2 (neutrino) part
  f[4]=(f2*c2-fL+xf3)/x2; // *y2 (anti-neutrino) part
}

int main()
{
  const double reps=.001;        // relative accuracy of the total Q2 integral calculation
  const double xeps=.0001;      // relative accuracy of the total X integral calculation
  //           =========
  const double GF=1.16637e-5;   // Fermi constant in GeV^-2
  const double GF2=GF*GF;       // Squared Fermi constant in GeV^-4
  const double MW=80.425;       // Mass of W-boson in GeV
  const double MW2=MW*MW;       // Squared mass of W-boson in GeV^2
  const double MW4=MW2*MW2;     // Quadro mass of W-boson in GeV^4
  const double hc2=38937932300.;// (hc)^2 in GeV^2*10^-38cm2 to convert GeV^-2 to 10^-38cm2
  const double pif=3.14159265*4;// 4pi
  const double sik=GF2*hc2/pif; // precalculated coefficient
  //const double mpi=.1349766;    // pi0 meson mass in GeV
  //const double mpi=.13957018;   // charged pi meson mass in GeV
  //const double mpi2=mpi*mpi;    // m_pi^2 in GeV^2
  //const double me=.00051099892; // electron mass in GeV
  //const double me2=me*me;       // m_e^2 in GeV^2
  //const double hme2=me2/2;      // .5*m_e^2 in GeV^2
  const double mmu=.105658369;  // mu meson mass in GeV
  const double mmu2=mmu*mmu;    // m_mu^2 in GeV^2
  const double hmmu2=mmu2/2;    // .5*m_mu^2 in GeV^2
  //const double mtau=1.777;      // tau meson mass in GeV
  //const double mtau2=mtau*mtau; // m_tau^2 in GeV^2
  //const double hmtau2=mtau2/2;  // .5*m_e^2 in GeV^2
  //const double mp=.93827203;    // proton mass in GeV
  //const double mn=.93956536;    // neutron mass in GeV
  //const double md=1.87561282;   // deuteron mass in GeV
  const double MN=.931494043;   // Nucleon mass (inside nucleus, atomic mass unit, GeV)
  //const double MN=(mn+mp)/2;    // Nucleon mass (mean free) in GeV
  //const double MD=1.232;        // proton mass in GeV
  //const double mp2=mp*mp;       // m_p^2 in GeV^2
  const double MN2=MN*MN;       // M_N^2 in GeV^2
  const double dMN=MN+MN;       // 2*M_N in GeV
  const double dMN2=MN2+MN2;    // 2*M_N^2 in GeV^2
  const double fMN2=dMN2+dMN2;  // 4*M_N^2 in GeV^2
  //const double EminE=me+me2/dMN;// Threshold for muon production
  const double EminMu=mmu+mmu2/dMN;  // Threshold for muon production
  //const double EminTau=mmu+mmu2/dMN; // Threshold for muon production
  //
  const double mc=.261;           // parameter of W2>(M_N^2+M_D^2)/2 cut for QuasiEl/Delta
  //const double mc=mpi;          // parameter of W>M+mc cut for Quasi-Elastic/Delta
  const double mcV=(dMN+mc)*mc; // constant of W>M+mc cut for Quasi-Elastic
  //std::ofstream fileNuMuX("NuMuXQ2.out", std::ios::out);
  //fileNuMuX.setf( std::ios::scientific, std::ios::floatfield );
  // _____ Begin of Test Area
  //Genfun::LogGamma logGamma;
  //double n=4.9;
		//double g=exp(logGamma(n));
  //G4cout<<"Gamma("<<n<<") = "<<g<<G4endl;
  // ^^^^^ End of Test Area
  //
  double f[5];                    // A working array
  int    A=12;                    // Neucleus for which calculations should be done
  double lEnuMin=0;               // LogLog of Minimum energy of neutrino
  double lEnuMax=std::log(1.+std::log(300./EminMu)); // LogLog of MaximumEnergy of neutrino
  int    nE=63;
  double dlE=(lEnuMax-lEnuMin)/nE;
  lEnuMin+=dlE/10;
  lEnuMax+=dlE/5;
  G4cout<<"Emin="<<EminMu<<",lEi="<<lEnuMin<<",lEa="<<lEnuMax<<",dlE="<<dlE<<G4endl;
  for(double lEnu=lEnuMin; lEnu<lEnuMax; lEnu+=dlE)
		{
    double Enu=std::exp(std::exp(lEnu)-1.)*EminMu; // Energy of neutrino/anti-neutrino
    double dEnu=Enu+Enu;          // doubled energy of nu/anu
    double Enu2=Enu*Enu;          // squared energy of nu/anu
    double Emu=Enu-mmu;           // Free Energy of neutrino/anti-neutrino
    double Emu2=Emu*Emu;          // squared energy of nu/anu
    double ME=Enu*MN;             // M*E
    double dME=ME+ME;             // 2*M*E
    double DIStsig=1.;            // Total curent DIS cross-section to be integrated
    double DISmsig=1.e20;         // Total remembered DIS cross-section
    double dEMN=(dEnu+MN)*ME;
    double MEm=ME-hmmu2;
    double sqE=Enu*std::sqrt(MEm*MEm-mmu2*MN2);
    double E2M=MN*Enu2-(Enu+MN)*hmmu2;
    double ymax=(E2M+sqE)/dEMN;
    double ymin=(E2M-sqE)/dEMN;
    double rmin=1.-ymin;
    double rhm2E=hmmu2/Enu2;
    double Q2min=(Enu2+Enu2)*(rmin-rhm2E-std::sqrt(rmin*rmin-rhm2E-rhm2E));
    double Q2max=dME*ymax;
    int    nQ2=8;
    //G4cout<<"*** E="<<Enu<<", Q2i="<<Q2min<<" < Q2a="<<Q2max<<", yi="<<ymin<<" < ya="
    //      <<ymax<<G4endl;
    while(std::fabs(DIStsig-DISmsig)/DIStsig>reps)
    {
      DISmsig=DIStsig;
      DIStsig=0.;
      nQ2*=2;
      double dQ2=(Q2max-Q2min)/nQ2;
      for(double Q2=Q2min+dQ2/2; Q2<Q2max; Q2+=dQ2)
				  {
        double DISxint=1.;        // Curent DIS x-integral
        double DISmint=1.e20;     // Remembered DIS x-integral
        double Q2M=Q2+MW2;
        double dik=MW4/Q2M/Q2M;
        double qmc=Q2+mcV;
        double lXQES=std::log((std::sqrt(qmc*qmc+Q2*fMN2)-qmc)/dMN2); // QuasielastBoundary
        //double lXQES=log(Q2/(Q2+mcV)); // Quasielastic boundary (W=MN+m_c)
        //double xN=Q2/dME;
        double xN=Q2/MN/(Emu+std::sqrt(Emu2+Q2));
        //double lXmin=log(xN/ymax);
        double lXmin=std::log(xN);
        // ****** QE ********
        if(lXQES>lXmin) lXmin=lXQES;   // A cut which leaves only QES
        // *** End of QE^^^^^
        double lXmax=0.;               // QES is in DIS
        //double lXmax=lXQES;          // Cut off quasielastic
        int    nX=8;
        while(std::fabs(DISxint-DISmint)/DISxint>xeps)
        {
          DISmint=DISxint;
          DISxint=0.;
          nX*=2;
          double dlX=(lXmax-lXmin)/nX;
          for(double lX=lXmin+dlX/2; lX<lXmax; lX+=dlX)
				      {
            getFun(A, lX, Q2, f);
            DISxint+=f[0]+f[0]+xN*(f[1]+f[1]+xN*f[3]); // neutrino
            //DISxint+=f[0]+f[0]+xN*(f[2]+f[2]+xN*f[4]); // anti-neutrino
            //G4cout<<f[0]<<","<<f[1]<<","<<f[2]<<","<<f[3]<<","<<f[4]<<G4endl;
          }
          DISxint*=dlX;
          //G4cout<<"--- E="<<Enu<<" --- Q2="<<Q2<<" --- nX="<<nX<<", iX="<<DISxint
          //      <<", mX="<<DISmint<<", rX="<<(DISxint-DISmint)/DISxint<<G4endl;
        }
        //G4cout<<"(E="<<Enu<<"), Q2="<<Q2<<", I="<<DISxint/dik/dik<<G4endl;
        DIStsig+=DISxint*dik;
      }
      DIStsig*=dQ2;
      //G4cout<<"=== E="<<Enu<<" ===> nQ="<<nQ2<<", iQ="<<DIStsig<<", mQ="<<DISmsig
						//      <<", rQ="<<(DIStsig-DISmsig)/DIStsig<<G4endl;
    }
    //===== tot/qe choice ====
    //DIStsig*=sik/Enu;
    //G4cout<<"***total*** E="<<Enu<<",sig/E= "<<DIStsig<<G4endl;
    //...................
    DIStsig*=sik;
    G4cout<<"***qelas*** E="<<Enu<<",sig= "<<DIStsig<<G4endl;
    //===== End of the choice
		} // End of the Enery LOOP
  // int np=0;
  //for(int m=0; m<2; m++)
  //{
  //  //fileNuMuX<<"  static const G4double SH"<<n<<"[nH]={"<<G4endl<<"    ";
  //  //G4cout<<"**** A_high="<<m<<G4endl;
  //  np=0;
  //  int nC=14;
	 //  for(G4int en=0; en<nC; en++)
	 //  {
  //    //G4double sig=1.;
  //    np++;
  //    //if(np==7)           // Write by 7 number in brackets
  //    //{
  //    //  if(en==nC-1) fileNuMuX<<sig<<"};"<<G4endl;
  //    //  else         fileNuMuX<<sig<<","<<G4endl<<"    ";
	 //    //}
  //    //else           fileNuMuX<<sig<<",";
  //    //if(np==7) np=0;
  //  } // End of the point LOOP
  //} // End of the isotop LOOP
  return EXIT_SUCCESS;
}
