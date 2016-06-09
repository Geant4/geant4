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

// Solution for randomization of the high order polinoms is tabulated

double poli(double x, int n) // f=n*x^(n-1)-(n-1)*x^n=R(0-1), n>2, as n=2 has analSolution
{
  double p=x;
  int n1=n-1;
  for(int i=1; i<n1; i++) p*=x;
  return p*(n-n1*x);
}

double dpoli(double x, int n) // df=n*(n-1)*x^(n-2)-(n-1)*n*x^(n-1), n>2 (n=2 has analSol)
{
  double p=x;
  int n1=n-1;
  if(n>3) for(int i=2; i<n1; i++) p*=x;
  return p*n*n1*(1.-x);
}
// One can define for tabulation more similar functions

int main()
{
  const double eps=.0000001;      // relative accuracy of the integral calculation
  //           =============
  int nSub=10;
  double dSub=1./nSub;
  for(int n=3; n<12; n++)
		{
    G4cout<<"n="<<n<<G4endl;
    G4double s=0.;
    for(int i=1; i<nSub; i++)
    {
      s+=dSub;
      G4double r=std::pow(s,(n-1.5)/1.5);
      G4double x=0.5;
		    if(n==3)
		    {
        if    (r==0.5) x=0.5;
        else if(r<0.5) x=std::sqrt(r+r)*(.5+.1579*(r-.5));
        else           x=1.-std::sqrt(2.-r-r)*(.5+.1579*(.5-r));
      }
      else
		    {
        G4int n1=n-1;
        G4double r1=n1;
        G4double r2=r1-1.;
        G4double rr=r2/r1;
        G4double rp=std::pow(rr,n1);
        G4double p2=rp+rp;
        if  (r==rr)  x=p2;
        else
				    {
          if(r<rr)
          {
								    G4double pr=0.;
								    G4double pra=0.;
            if(n>7)
								    {
              if(n>9)
								      {
                if(n>10)                         // >10(11)
                {
                  pr=.614/std::pow((n+1+1.25),.75);
                  pra=.915/std::pow((n+1+6.7),1.75);
                }
												    else                             // 10
                {
                  pr=.09945;
                  pra=.00667;
                }
              }
              else
								      {
                if(n>8)                          // 9
                {
                  pr=.1064;
                  pra=.00741;
                }
												    else                             // 8
                {
                  pr=.11425;
                  pra=.00828;
                }
              }
            }
            else
								    {
              if(n>5)
								      {
                if(n>6)                          // 7
                {
                  pr=.12347;
                  pra=.00926;
                }
												    else                             // 6
                {
                  pr=.13405;
                  pra=.01027;
                }
              }
              else
								      {
                if(n>4)                          // 5
                {
                  pr=.1454;
                  pra=.01112;
                }
												    else                             // 4
                {
                  pr=.15765;
                  pra=.00965;
                }
              }
            }
            x=std::pow((r/p2),(1.-rr+pra))*(rr+pr*(r-p2));
          }
          else
          {
								    G4double sr=0.;
            if(n>7)
								    {
              if(n>9)
								      {
                if(n>10) sr=.86/(n+1+1.05);      // >10(11)
												    else     sr=.0774;               // 10
              }
              else
								      {
                if(n>8) sr=.0849;                // 9
												    else    sr=.0938;                // 8
              }
            }
            else
								    {
              if(n>5)
								      {
                if(n>6) sr=.1047;                // 7
												    else    sr=.1179;                // 6
              }
              else
								      {
                if(n>7) sr=.1339;                // 5
												    else    sr=.15135;               // 4
              }
            }
            x=1.-std::sqrt((1.-r)/(1.-p2))*(1.-rr+sr*(p2-r));
          }
        }
      }
      G4double dx=x;
      G4double f=poli(x,n)-r;
      G4int it=0;
      while(std::fabs(f/r)>eps)
						{
        it++;
        G4double df=dpoli(x,n);
        //G4cout<<"n="<<n<<", r="<<r<<", f="<<f<<", d="<<df<<", x="<<x<<", i="<<it<<G4endl;
        x-=f/df;
        f=poli(x,n)-r;
      }
      dx=std::fabs(dx-x);
      G4double d=std::fabs(f/r);
      G4cout<<"n="<<n<<",r="<<r<<", Final: f="<<d<<",x="<<x<<",d="<<dx<<",it="<<it<<G4endl;
    } // End of loop over r
  } // End of loop over n
  return 0;
}
