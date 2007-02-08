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
#include "Randomize.hh"
#include "G4Timer.hh"
#include "time.h"

// All calculations have been done for C12 nucleus

double anuX(double E, double r, double p) // (E is in GeV, r=Q2/Q2max, p=1+anuE(E))
{
		double	E2=E*E;
  double y=p-r;
  double p3=(13.88+.9373*(1.+.000033*E2)*std::sqrt(E))/(1.+(10.12+1.532/E2)/E);
  return std::pow(y,6)*(r+p3)/(p3*r+1.);
}

int main()
{
  G4double PI=3.14159265;        // \pi
  G4double A=.5;                 // half side of the tested cube
  G4double R=.5;                 // radius of the tested sphere
  G4double B=1.;                 // half width of the mother cube
  G4double S=6*(B+B)*(B+B);      // 6 sides 2x2
  G4double SC=6*(A+A)*(A+A);     // The target surface of the test cube
  G4double SS=4*PI*R*R;          // The target surface of the test sphere
  G4double c=0.;                 // summ for the cube							
  G4double s=0.;                 // summ for the sphere					
  G4double m=0.;                 // summ for the mother cube
  G4double cf=0.;                // reversed length summ for the cube							
  G4double sf=0.;                // reversed length summ for the sphere					
  G4double mf=0.;                // reversed length summ for the mother cube
  G4double cr=0.;                // length summ for the cube							
  G4double sr=0.;                // length summ for the sphere					
  G4double mr=0.;                // length summ for the mother cube
  G4double cm=B;                 // min length for the cube							
  G4double sm=B;                 // min length for the sphere					
  G4double mm=B;                 // min length for the mother cube
  G4int nEv=10000000;            // calculation statistics
  G4double f=.99;                // scale factor
  G4double Af=A*f;               // Scaled A
  G4double Rf=R*f;               // Scaled R
  G4double Bf=B*f;               // Scaled B
  G4Timer* timer = new G4Timer();
  timer->Start();
  for(G4int i=0; i<nEv; i++)
		{
    // Randomize coordinates within (+/-1,+/-1,+/-1) cube (S=6+4=24)
    G4double x=G4UniformRand();
    x=B*(1.-x-x);                    // Can be x=B*x, but B=1 is a hilf width of mother cube
    G4double y=G4UniformRand();
    y=B*(1.-y-y);
    G4double z=G4UniformRand();
    z=B*(1.-z-z);
    G4double r2=x*x+y*y+z*z;
    G4double r=std::sqrt(r2);
    G4double w=1./r;
    G4double ax=std::fabs(x);
    G4double ay=std::fabs(y);
    G4double az=std::fabs(z);
    if(r<R && r>Rf)
    {
      s+=1.;
      sf+=w;
      sr+=r;
      if(r<sm) sm=r;
    }
    if(ax<A && ay<A && az<A && (ax>Af || ay>Af || az>Af))
    {
      c+=1.;
      cf+=w;
      cr+=r;
      if(r<cm) cm=r;
    }
    if(ax>Bf || ay>Bf || az>Bf)
    {
      m+=1.;
      mf+=w;
      mr+=r;
      if(r<mm) mm=r;
    }
  }
  G4double N=S*mm/m;
  //G4double N=S*m*m/mf;
  //G4double N=S/mf;
  timer->Stop();
  G4cout<<"SurfaceCalc: nE="<<nEv<<", time="<<timer->GetUserElapsed()<<" s, f="<<f<<G4endl;
  delete timer;
  G4cout<<"SurfaceCalculation: Cube="<<N*c/cm/SC<<", Sphere="<<N*s/sm/SS<<G4endl;
  //G4cout<<"SurfaceCalculation: Cube="<<N*cf/SC<<", Sphere="<<N*sf/SS<<G4endl;
  //G4cout<<"SurfaceCalculation: Cube="<<N*cr/c/c/SC<<", Sphere="<<N*sr/s/s/SS<<G4endl;
}
