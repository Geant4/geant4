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

double nuE(double E) { return .9673/(1.+.323/E/E)/std::pow(E,.78);} // (E is in GeV)

double nuX(double E, double r, double p) // (E is in GeV, r=Q2/Q2max, p=1.+nuE(E))
{
  double y=p-r;
  double p3=(.3088+.0012352*E)/(1.+1.836/E/E);
  return std::pow(y,6)*(r+p3)/(p3*r+1.);
}

double anuE(double E) { return .875/(1.+.2977/E/E)/std::pow(E,.78);} // (E is in GeV)

double anuX(double E, double r, double p) // (E is in GeV, r=Q2/Q2max, p=1+anuE(E))
{
		double	E2=E*E;
  double y=p-r;
  double p3=(13.88+.9373*(1.+.000033*E2)*std::sqrt(E))/(1.+(10.12+1.532/E2)/E);
  return std::pow(y,6)*(r+p3)/(p3*r+1.);
}

int main()
{
  const double eps=.00000001;    // relative accuracy of the integral calculation
  const double mmu=.105658369;   // mu meson mass in GeV
  const double mmu2=mmu*mmu;     // m_mu^2 in GeV^2
  const double MN=.931494043;    // Nucleon mass (inside nucleus, atomic mass unit, GeV)
  const double Emin=mmu+mmu2/(MN+MN); // the threshold energy in GeV
  const double Emax=390.;        // the maximum energy in GeV
  const double lEmi=std::log(Emin);   // logarithm of the threshold energy in GeV
  const double lEma=std::log(Emax);   // logarithm of the maximum energy in GeV
  const int    power=7;          // power for the magic variable (E-dependent)
  const double pconv=1./power;   // conversion power for the magic variable
  //           =========
  const int niX=20;              // number of points for X table
  const int liX=niX-1;           // the last index for X table
  const int niE=20;              // number of points for E table
  double Xl[niE][niX];           // reversed table
  double inl[niE][niX];          // direct table
  // ----------------- nu/anu switch ------------------------------
  //bool nu=true;
  bool nu=false;
  // --------------------------------------------------------------
  double dE=(lEma-lEmi)/niE;
  double hE=dE/2;
  int ne=0;
  for(double len=lEmi+hE; len<lEma; len+=dE)
  {
    //G4cout<<"log(E)="<<len<<G4endl;
    double en=std::exp(len);
    double shift=0.;
    if(nu) shift=1.+nuE(en);
    else   shift=1.+anuE(en);
    double Xma=std::pow(shift,power);
    double Xmi=std::pow((shift-1.),power);
    int    nX=8;
    double DISmsig=0.;
    double DIStsig=1.;
    while(std::fabs(DIStsig-DISmsig)/DIStsig>eps)
    {
      DISmsig=DIStsig;
      DIStsig=0.;
      nX*=2;
      double dX=(Xma-Xmi)/nX;
      double hX=dX/2;
      for(double X=Xmi+hX; X<Xma; X+=dX)
				  {
        double r=shift-std::pow(X,pconv); // the same for nu and anu
        if(nu) DIStsig+=nuX(en,r,shift);  // neutrino
        else   DIStsig+=anuX(en,r,shift); // anti-neutrino
      }
      DIStsig*=dX;
      //G4cout<<"E="<<en<<",nX="<<nX<<",i="<<DIStsig<<",m="<<DISmsig<<",r="
      //      <<(DIStsig-DISmsig)/DIStsig<<G4endl;
		  }
    //G4cout<<"E="<<en<<", Total: nX="<<nX<<", Int="<<DIStsig<<G4endl;
    // ***************** Calculate the reversed table *****************
    double dInt=DIStsig/niX;             // Step for the integral
    DIStsig=0.;
    DISmsig=0.;
    nX*=2;
    double dX=(Xma-Xmi)/nX;      // New step for the multiplied number of nodes
    dInt/=dX;                    // To avoid multiplication (nX is fixed)
    double sum=dInt;
    int nn=0;
    double hX=dX/2;
    G4cout<<"E="<<en<<", Xl_min="<<Xmi<<G4endl;
    for(double X=Xmi+hX; X<Xma; X+=dX)
		  {
      double r=shift-std::pow(X,pconv); // the same for nu and anu
      if(nu) DIStsig+=nuX(en,r,shift);  // neutrino
      else   DIStsig+=anuX(en,r,shift); // anti-neutrino
      if(DIStsig>sum+eps)
		  		{
        inl[ne][nn]=sum*dX;
        Xl[ne][nn]=X-(DIStsig-sum)*dX/(DIStsig-DISmsig);
        G4cout<<"E="<<en<<", sum="<<inl[ne][nn]<<", Xl["<<nn<<"]="<<Xl[ne][nn]<<G4endl;
        nn++;
        sum+=dInt;
      }
      DISmsig=DIStsig;
    }
    inl[ne][nn]=sum*dX;
    Xl[ne][nn]=Xma;
    G4cout<<"E="<<en<<", sum="<<inl[ne][nn]<<", Xl["<<nn<<"]="<<Xma<<G4endl;
    // ************* The following is just a test (better than .5%) ************
    //for(int i=1; i<niX; i++)
    //{
    //  DIStsig=0.;
    //  double Xm=Xl[ne][i];
    //  double dX=(Xm-Xmi)/nX;
    //  double hX=dX/2;
    //  for(double X=Xmi+hX; X<Xm; X+=dX)
		  //		{
    //    double r=shift-pow(X,pconv); // the same for nu and anu
    //    if(nu) DIStsig+=nuX(en,r,shift);  // neutrino
    //    else   DIStsig+=anuX(en,r,shift); // anti-neutrino
    //  }
    //  DIStsig*=dX;
    //  G4cout<<"i="<<i<<", v="<<DIStsig<<", d="<<fabs(DIStsig-inl[ne][i])/DIStsig<<G4endl;
    //}
    // ***************** Calculate the direct table *****************
    dX=(Xma-Xmi)/niX;
    double nor=inl[ne][liX]/niX;
    //G4cout<<"E="<<en<<", i=0, Xm="<<Xmi<<", I=0."<<G4endl;
    for(int i=1; i<=niX; i++)
    {
      DIStsig=0.;
      double Xm=Xmi+dX*i;
      double rX=(Xm-Xmi)/nX;
      double hX=rX/2;
      for(double X=Xmi+hX; X<Xm; X+=rX)
		  		{
        double r=shift-std::pow(X,pconv); // the same for nu and anu
        if(nu) DIStsig+=nuX(en,r,shift);  // neutrino
        else   DIStsig+=anuX(en,r,shift); // anti-neutrino
      }
      DIStsig*=rX/nor;
      //////////G4cout<<"E="<<en<<", i="<<i<<", Xm="<<Xm<<", I="<<DIStsig<<G4endl;
    }
    ne++;
  } // End of the big loop over log(E)
  G4cout<<"End"<<G4endl;
  return 0;
}
