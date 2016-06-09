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

// Calculation of electron spectra integrals for randomization of bound mu->e+nu+nu decay

double espec(double Z , double T) // T=E_kin in MeV
{
  static double mZ=0.,mp,mc,mep,mem,me,mpp;
  if(Z!=mZ)
  {
    double Z2=Z*Z;
    double Z4=Z2*Z2;
    double az=std::log(Z);
    mZ=Z;
    mp=(11.65+330./Z)/(1.+7.6e-17*Z4*Z4);
    mc=52.83-.2257*az*az;
    mep=.674;
    mem=.1431+Z2/15080.;
    me=5.54*std::exp(-.0327*std::pow(Z,1.241));
    mpp=2.253*std::exp(-Z2/8922.7);
  		G4cout<<"z="<<Z<<":p="<<mp<<",c="<<mc<<",ep="<<mep<<",em="<<mem<<",e="<<me<<",pp="<<mpp
          <<G4endl;
		}
  double p=mp;
  double c=mc;
  double ep=mep;
  double em=mem;
  double e=me;
  double pp=mpp;
  double res=std::pow(T,pp)/(1.+e*std::exp(em*std::pow(T,ep)))/(1.+std::pow(T/c,p));
  //G4cout<<"T="<<T<<", R="<<.0001389*res<<G4endl;
  return res;
}
int main()
{
  const int niX=201;             // number of points for randomization table
  const int liX=niX-1;           // the last index for randomization table (?)
  //const int niE=1000000;         // number of points for integration table
  const int niE=1000;           // number of points for integration table
  // **********************************************************************
  double Z=20.;                  // Z of the Element **********************
  // **********************************************************************
  double X[niX];                 // randomization table
  double in[niE];                // integration table
  double Emax=80.;               // Max kin E in MeV
  double dE=Emax/niE;            // Step in the integration table
  double spP=0.;                 // Previous spectrum value [sp(0)=0]
  in[0]=0.;                      // Not yet normed collected integrals
  double e=0.;                   // running value of energy
  double inv=0.;                 // Integral cumulative value
  for(int i=1; i<niE; i++)
  {
    e+=dE;
    double sp=espec(Z,e);
    inv+=spP+sp;                 // integrate with the Gauss method
    in[i]=inv;
  		//G4cout<<"sp="<<sp<<",in[="<<i<<"]="<<inv<<G4endl;
    spP=sp;
  } // End of the big loop over log(E)
  for(int j=1; j<niE; j++) in[j]/=inv; // Normalize the integrals by a unit
  //for(int n=1; n<niE; n++) G4cout<<"1-in["<<n<<"]="<<1.-in[n]<<G4endl;
  double dx=1./liX;              // The last in the table is MAXE by definition
  double d=dx;                   // Start searching with 1, because x(0)=0
  int m=1;                       // index in the randomization table
  double oi=0.;                  // previouse integral value
  X[0]=0.;                       // Interpolated randomization table (energy)
  e=0.;
  for(int k=1; k<niE; k++)
		{
    e+=dE;
    double ci=in[k];
    if(ci>=d)
				{
      if(ci==d) X[m]=e;
      else X[m]=e-(ci-d)*dE/(ci-oi);
      m++;
      d+=dx;                     // start searching for another point
    }
				oi=ci;
  }
  for(int k=0; k<liX; k++)
		{
    G4cout<<X[k]<<",";
    if(!((k+1)%10))G4cout<<G4endl;
  }
  return 0;
}
