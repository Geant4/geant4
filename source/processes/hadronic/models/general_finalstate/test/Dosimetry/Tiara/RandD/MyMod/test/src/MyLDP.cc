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
// $Id: MyLDP.cc,v 1.1 2003-10-08 12:32:19 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "MyLDP.hh"


const G4double MyLDP::ConstEvapLevelDensityParameter = 0.125*(1./MeV);
const G4double MyLDP::alpha = 0.076*(1./MeV);
const G4double MyLDP::beta = 0.180*(1./MeV);
const G4double MyLDP::beta1 = 0.157*(1./MeV);
const G4double MyLDP::gamma = 0.47*(1./MeV);
const G4double MyLDP::Bs = 1.0;


MyLDP::MyLDP(const MyLDP &right)
{
    G4Exception("G4EvaporationLevelDensityParameter::copy_constructor meant to not be accessable");
}


const MyLDP & MyLDP::operator=(const MyLDP &right)
{
    G4Exception("G4EvaporationLevelDensityParameter::operator= meant to not be accessable");
    return *this;
}


G4bool MyLDP::operator==(const MyLDP &right) const
{
    return false;
}

G4bool MyLDP::operator!=(const MyLDP &right) const
{
    return true;
}

G4double MyLDP::LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const 
{
    G4int N = A - Z;

    // Asymptotic Level Density Parameter
    G4double AsymptoticLDP = alpha*G4double(A) + beta*pow(G4double(A),2./3.) + beta1*pow(G4double(A),1./3.);
	
    // Shape of the LDP U dependence
    G4double exponent = -gamma*U/pow(G4double(A),1./3.);
    G4double f = 1.;
    if (exponent > -300.) f -= exp(exponent);
	
    // Level Density Parameter
    G4double a = AsymptoticLDP*(1. + ShellCorrection(Z,N)*f/U);
    return a;
    if(U>50*MeV) return a;
    G4double aNew,a1=a,T;
    for(unsigned i=0;i<100;i++){
      T = sqrt(U/a);
      //T = T*T*pow(G4double(A),2./3.)/441.;
      exponent = T*T*pow(G4double(A),2./3)/441;
      f = 1;
      if(exponent < 300) f = exp(-exponent);
      aNew = pow(0.7143*(1.-0.4*f),1.7);
      if(fabs(aNew-a)/aNew < 0.005){
	a = aNew;
	break;
      }
      a = aNew;
    }
    G4cout<<a<<G4endl;
    return a;
}

