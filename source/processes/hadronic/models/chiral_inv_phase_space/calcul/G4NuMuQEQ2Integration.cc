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
// G4 Tools program: NuMu DIS(Q2) fixed step integration
// .....................................................
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-2005
// 
//=====================================================================
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include "G4ios.hh"
//#include <CLHEP/GenericFunctions/LogGamma.hh>

// All calculations have been done for C12 nucleus

double nuQ2(double Q2) // (Q2 is in GeV^2)
{
  //  x      *     *     *     *     *     *     *     *      *     *
  // -1-    -2-   -3-   -4-   -5-   -6-   -7-   -8-   -9-   -10-  -11-
  // 2.068 2.634 4.333 3.000 10.59 .0011 18.74 1.484 134950 .0755 4.5
  G4double y=Q2/3;
  G4double z=std::pow((Q2/2.634),4.333);
		// The first (1.+Q2)**4.5 fuctor is just normalization, which payed back by 1/(1+Q2)^3.5
  G4double f=std::pow((1.+Q2),4.5)*
				(1.+z*(1.+18.74*std::exp(-Q2/1.484)-134950*std::exp(-Q2/.0755)))/std::pow((1.+y+.0011*y*y),10.59);
  //G4cout<<"Q2="<<Q2<<", y="<<y<<", z="<<z<<", f="<<f<<G4endl;
  return f;
}

double anuQ2(double Q2)
{
  //    x     *     *     *     *     *     *     *     -     *
  //   -1-   -2-   -3-   -4-   -5-   -6-   -7-   -8-   -9-   -10-
  // .2000 .0100 .3000 .4000 3.900 .0150 .2000 6.500 0.000 .00000001
  G4double  y=Q2/.4;
  G4double  y2=y*y;
  G4double  z=std::pow((Q2/.01),-.3);
		// The first (1.+Q2)**6.5 fuctor is just normalization, which payed back by 1/(1+Q2)^5.5
  G4double f=(1.+z)*std::pow((1+Q2),6.5)*std::pow(Q2,-.2)/std::pow((1.+y+(.015+.00000001*y2)*y2),3.9);
  //G4cout<<"Q2="<<Q2<<", y="<<y<<", z="<<z<<", f="<<f<<G4endl;
  return f;
}

int main()
{
  const double eps=.000001;      // relative accuracy of the integral calculation
  //           =========
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
  //const double dMN2=MN2+MN2;    // 2*M_N^2 in GeV^2
  //const double fMN2=dMN2+dMN2;  // 4*M_N^2 in GeV^2
  //const double EminE=me+me2/dMN;// Threshold for muon production
  const double Emin=mmu+mmu2/dMN;  // Threshold for muon production
  //const double EminTau=mmu+mmu2/dMN; // Threshold for muon production
  //
  //const double mc=.3;           // parameter of W>M+mc cut for Quasi-Elastic/Delta
  //const double mc=mpi;          // parameter of W>M+mc cut for Quasi-Elastic/Delta
  //const double mcV=(dMN+mc)*mc; // constant of W>M+mc cut for Quasi-Elastic
  //std::ofstream fileNuMuX("NuMuXQ2.out", std::ios::out);
  //fileNuMuX.setf( std::ios::scientific, std::ios::floatfield );
  // _____ Begin of Test Area
  //Genfun::LogGamma logGamma;
  //double n=4.9;
		//double g=exp(logGamma(n));
  //G4cout<<"Gamma("<<n<<") = "<<g<<G4endl;
  // ^^^^^ End of Test Area
  //
  const int niQ2=100;           // Number of points in the Q2 integration
  const int diQ2=niQ2+1;        // Dimention for arrays
  double Xl[diQ2];
  double inl[diQ2];
  double Enu=Emin;              // Initial E=minE  
  double Emo=Emin+Emin;         // First new
  double Emv=0.;
  double Q2i=1.;
  double Q2o=1.;
  double Q2v=0.;
		int nit=0;
  while (std::fabs(Q2v-Q2o)>eps && nit<10)
  {
    nit++;
    double Emm=Emo;               // Candidate for the next point is predefined
    double dEnu=Enu+Enu;          // doubled energy of nu/anu
    double Enu2=Enu*Enu;          // squared energy of nu/anu
    double ME=Enu*MN;             // M*E
    double dEMN=(dEnu+MN)*ME;
    double MEm=ME-hmmu2;
    double sqE=Enu*std::sqrt(MEm*MEm-mmu2*MN2);
    double E2M=MN*Enu2-(Enu+MN)*hmmu2;
    //double ymax=(E2M+sqE)/dEMN;
    double ymin=(E2M-sqE)/dEMN;
    double rmin=1.-ymin;
    double rhm2E=hmmu2/Enu2;
    Q2i=(Enu2+Enu2)*(rmin-rhm2E-std::sqrt(rmin*rmin-rhm2E-rhm2E)); // minQ2 for Enu
    //G4cout<<"MinSearch: E="<<Enu<<": Q2="<<Q2i<<G4endl;
    if(Emv>0.0001) // Not initialization (not first two) Use three points
				{
      //double d=Enu*Emo*(Enu-Emo)+Enu*Emv*(Emv-Enu)+Emo*Emv*(Emo-Emv);
      double a=Q2i*Emo-Enu*Q2o+Q2v*Enu-Emv*Q2i+Q2o*Emv-Emo*Q2v;
						double b=Enu*Enu*(Q2o-Q2v)+Emo*Emo*(Q2v-Q2i)+Emv*Emv*(Q2i-Q2o);
      //double a
				  Emm=-b/(a+a);
      if(Q2v<Q2i)                 // swap to make Q2v the biggest
				  {
        double q2=Q2v;
        double en=Emv;
        Q2v=Q2i;
        Emv=Enu;
        Q2i=q2;
        Enu=en;
      }
      if(Q2v<Q2o)                 // swap to make Q2v the biggest
				  {
        double q2=Q2v;
        double en=Emv;
        Q2v=Q2o;
        Emv=Emo;
        Q2o=q2;
        Emo=en;
      }
    }
    else if(Q2o>0.0000000001) Emm=Enu+Emin; // the second step
    //G4cout<<"***Enu="<<Enu<<", Emm="<<Emm<<", Emin="<<Emin<<G4endl;
    Emv=Emo;
    Q2v=Q2o;
    Emo=Enu;
    Q2o=Q2i;
    Enu=Emm;
    //G4cout<<"___Q2i="<<Q2i<<": Q2o="<<Q2o<<", Q2v="<<Q2v<<G4endl;
  }
  double Q2min=Q2i-eps;
  double Q2max=600.;            // covers the calculated region
  // ----------------- nu/anu switch ------------------------------
  //bool nu=true;
  bool nu=false;
  // -----------------------------**************************-------
  // *************** Convert to the magic variable ****************
  double Xmin=0.;
  double Xmax=0.;
		if(nu)
  {
    Xmin=std::pow(1+Q2max,-3.5);
    Xmax=std::pow(1+Q2min,-3.5);
		}  
  else
  {
    Xmin=std::pow(1+Q2max,-5.5);
    Xmax=std::pow(1+Q2min,-5.5);
		}
  G4cout<<"Ei="<<Enu<<":Q2i="<<Q2min<<",Q2a="<<Q2max<<",Xmi="<<Xmin<<",Xma="<<Xmax<<G4endl;
  int    nQ2=8;                 // nitial #Of points for the overall integration
  double DISmsig=0.;
  double DIStsig=1.;
		double Q2=0.;                 // Prototype for the integration 
  while(std::fabs(DIStsig-DISmsig)/DIStsig>eps)
  {
    DISmsig=DIStsig;
    DIStsig=0.;
    nQ2*=2;
    double dX=(Xmax-Xmin)/nQ2;
    double hX=dX/2;
    for(double X=Xmin+hX; X<Xmax; X+=dX)
				{
      if(nu) // neutrino
						{
        Q2=std::pow(X,-.2857142); // -1./3.5
        DIStsig+=nuQ2(Q2);
      }
      else  // anti-neutrino
						{
        Q2=std::pow(X,-.1818182); // -1./5.5
        DIStsig+=anuQ2(Q2);
      }
    }
    DIStsig*=dX;
    //G4cout<<"n="<<nQ2<<",i="<<DIStsig<<",m="<<DISmsig<<",r="<<(DIStsig-DISmsig)/DIStsig
				//		    <<G4endl;
		}
  G4cout<<"Total: nQ2="<<nQ2<<", Int="<<DIStsig<<G4endl;
  // ***************** Calculate the reversed table *****************
  double dInt=DIStsig/niQ2;             // Step for the integral
  DIStsig=0.;
  DISmsig=0.;
  nQ2*=2;
  double dX=(Xmax-Xmin)/nQ2;
  dInt/=dX;                     // To avoid multiplication (nQ2 is fixed)
  double sum=dInt;
  int nn=1;
  Xl[0]=0.;
  double hX=dX/2;
  for(double X=Xmin+hX; X<Xmax; X+=dX)
		{
    if(nu) // neutrino
				{
      Q2=std::pow(X,-.2857142); // -1./3.5
      DIStsig+=nuQ2(Q2);
    }
    else  // anti-neutrino
				{
      Q2=std::pow(X,-.1818182); // -1./5.5
      DIStsig+=anuQ2(Q2);
    }
    if(DIStsig>sum+eps)
				{
      inl[nn]=sum*dX;
      Xl[nn]=X-(DIStsig-sum)*dX/(DIStsig-DISmsig);
      G4cout<<"sum="<<inl[nn]<<", Xl["<<nn<<"]="<<Xl[nn]<<G4endl;
      nn++;
      sum+=dInt;
    }
    DISmsig=DIStsig;
  }
  inl[nn]=sum*dX;
  Xl[nn]=Xmax;
  G4cout<<"sum="<<inl[nn]<<", Xl["<<nn<<"]="<<Xmax<<G4endl;
  // ************* The following is just a test (better than .5%) ************
  for(int i=1; i<=niQ2; i++)
  {
    DIStsig=0.;
    double Xm=Xl[i];
    double dX=(Xm-Xmin)/nQ2;
    double hX=dX/2;
    for(double X=Xmin+hX; X<Xm; X+=dX)
				{
      if(nu) // neutrino
						{
        Q2=std::pow(X,-.2857142); // -1./3.5
        DIStsig+=nuQ2(Q2);
      }
      else  // anti-neutrino
						{
        Q2=std::pow(X,-.1818182); // -1./5.5
        DIStsig+=anuQ2(Q2);
      }
    }
    DIStsig*=dX;
    G4cout<<"i="<<i<<", v="<<DIStsig<<", d="<<std::fabs(DIStsig-inl[i])/DIStsig<<G4endl;
  }
  // ***************** Calculate the direct table *****************
  dX=(Xmax-Xmin)/niQ2;
  double nor=inl[niQ2]/niQ2;
  for(int i=1; i<=niQ2; i++)
  {
    DIStsig=0.;
    double Xm=Xmin+dX*i;
    double rX=(Xm-Xmin)/nQ2;
    double hX=rX/2;
    for(double X=Xmin+hX; X<Xm; X+=rX)
				{
      if(nu) // neutrino
						{
        Q2=std::pow(X,-.2857142); // -1./3.5
        DIStsig+=nuQ2(Q2);
      }
      else  // anti-neutrino
						{
        Q2=std::pow(X,-.1818182); // -1./5.5
        DIStsig+=anuQ2(Q2);
      }
    }
    DIStsig*=rX/nor;
    G4cout<<"i="<<i<<", Xm="<<Xm<<", I="<<DIStsig<<G4endl;
  }
  G4cout<<"End"<<G4endl;
  DIStsig=0.;
  return 0;
}
