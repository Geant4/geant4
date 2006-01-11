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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: AnaEx02.cc,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
//      GEANT 4 - AnaEx01
//
// --------------------------------------------------------------------
//   Example of G4 fit usage.
//   See the README file within the same directory to have more infos.
// --------------------------------------------------------------------

#include "Randomize.hh"
#include "G4QFunFit.hh"
#include "AnaEx02Function.hh"

int main(int,char**)
{
  const G4int Npr= 6;                        // #of parameters
  G4int Nar= 4;                              // #of arguments
  G4int Nb = 3;                              // basis for the power
  G4int Npt= Nb;                             // #of points start
  for(G4int I=1; I<Nar; I++) Npt*=Nb;        // #of points result
  G4int pC = 0;                              // counter of points
  G4double B0 = 1.;                      // Base for X0 
  G4double B1 = 1.;                          // Base for X1 
  G4double B2 = 1.;                          // Base for X2 
  G4double B3 = 1.;                          // Base for X3 
  G4double D0 = 1.;                      // Step for X0 
  G4double D1 = 1.;                          // Step for X1 
  G4double D2 = 1.;                          // Step for X2 
  G4double D3 = 1.;                          // Step for X3 
  G4double AV[Npr] = {1.,2.,3.,4.,5.,6.};    // values of parameters
  G4double MI[Npr] = {.1,.1,.1,.1,.1,.1};    // min values of parameters
  G4double MA[Npr] = {9.,9.,9.,9.,9.,9.};    // max values of parameters
  std::vector<G4double>* X  = new std::vector<G4double>; // Arguments
  std::vector<G4double>* A  = new std::vector<G4double>; // Parameters
  std::vector<G4double>* ER = new std::vector<G4double>; // Errors(ranges) of parameters
  std::vector<G4double>* AI = new std::vector<G4double>; // Minimum A
  std::vector<G4double>* AA = new std::vector<G4double>; // Maximum A
  std::vector<G4double>* EX = new std::vector<G4double>; // Ex.data Npt*(Y,dY,X0,..,XN)
  G4double DY=0.4;
  for(G4int I=0; I<Npr; I++)
		{
    A->push_back(AV[I]);
    ER->push_back(.5);
    AI->push_back(MI[I]);
    AA->push_back(MA[I]);
  }
  for(G4int J=0; J<Nar; J++) X->push_back(0.); // activate arguments
  AnaEx02Function* fun = new AnaEx02Function(A, X, AI, AA);  
  // Fill Ex. data with arguments
  for(G4int I0=0; I0<Nb; I0++)
		{
    G4double C0=B0+D0*I0;
    fun->SetX(0, B0+D0*I0);
    for(G4int I1=0; I1<Nb; I1++)
		  {
      G4double C1=B1+D1*I1;
      fun->SetX(1, B1+D1*I1);
      for(G4int I2=0; I2<Nb; I2++)
		    {
        G4double C2=B2+D2*I2;
        fun->SetX(2, B2+D2*I2);
        for(G4int I3=0; I3<Nb; I3++)
		      {
          G4double C3=B3+D3*I3;
          fun->SetX(3, C3);
          G4double Y0=fun->Y();       // Without parameters -> only Y
          G4double dY=G4UniformRand();
          Y0+=dY+dY-1.;
          EX->push_back(Y0);          // Y
          EX->push_back(DY);          // dY
          EX->push_back(C0);          // X0 
          EX->push_back(C1);          // X1 
          EX->push_back(C2);          // X2 
          EX->push_back(C3);          // X3 
          pC++;
       }
      }
    }
  }
  G4cout<<"AnaEx02: #of points counter="<<pC<<" =? Npt="<<Npt<<G4endl;
  G4cout<<" A:";
  for(G4int I=0; I<Npr; I++)
		{
    if(!I) (*A)[I]+=.3*(G4UniformRand()-.5);
    else   (*A)[I]+=.1*(G4UniformRand()-.5);
    G4cout<<" "<<(*A)[I];
  }
  G4cout<<G4endl;
  G4QFunFit* funfit = new G4QFunFit(0,A,ER,EX,fun,Npt,AI,AA); // 0=Hi2 minimization
  G4double S= funfit->FIT(2,1,200,.01,10);
  //G4int res = funfit->GetENDFLG();
  G4int def = Npt-Npr;
  //G4cout<<"AnaEx02: Fit is made with MC="<<res<<", Hi2/ns="<<(S+S)/def<<G4endl;
  G4cout<<"AnaEx02: Fit is made with  Hi2/ns="<<(S+S)/def<<G4endl;
  std::vector<G4double>* SI=funfit->GetSIG();  
  for(G4int I=0; I<Npr; I++)
		{
    G4cout<<I<<", A="<<(*A)[I];
    if((*ER)[I])G4cout<<", D="<<(*ER)[I];
    else G4cout<<", IS FIXED ON BOUNDARY";
    G4cout<<", SIG="<<(*SI)[I]<<G4endl;
  }
  delete funfit;
  delete fun;
  return 0;
}
