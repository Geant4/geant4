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
// $Id: PreCompoundParameters.cc,v 1.1 2003-10-08 12:32:13 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#include "PreCompoundParameters.hh"


PreCompoundParameters PreCompoundParameters::thePreCompoundParameters;

PreCompoundParameters * PreCompoundParameters::GetAddress()
{ return &thePreCompoundParameters; }

#define alpha 0.076/MeV
#define beta  0.180/MeV
#define beta1  0.157/MeV
#define gamma  0.47/MeV
#define Bs  1.0
#include "G4CameronTruranHilfShellCorrections.hh"
G4double PreCompoundParameters::GetLevelDensity(const G4Fragment& fragm)
{
  G4int A = fragm.GetA(),Z = fragm.GetZ(),N = A - Z;
  G4double U = fragm.GetExcitationEnergy();

    // Asymptotic Level Density Parameter
  G4double AsymptoticLDP = alpha*G4double(A) + beta*pow(G4double(A),2./3.) + beta1*pow(G4double(A),1./3.);
  // Shape of the LDP U dependence
  G4double exponent = -gamma*U/pow(G4double(A),1./3.);
  G4double f = 1.;
  if (exponent > -300.) f -= exp(exponent);
        
    // Level Density Parameter
  G4double a = AsymptoticLDP*(1. + G4CameronTruranHilfShellCorrections::GetInstance()->GetShellCorrection(Z,N)*f/U);
  if(U>50) return a;
  G4double aNew,a1=a,T;
  for(unsigned i=0;i<100;i++){
    T = sqrt(U/a);
    //T = T*T*pow(G4double(A),2./3.)/441.;
    exponent = T*T*pow(A,2./3.)/441;
    f = 1;
    if(T < 300) f = exp(-exponent);
    aNew = a1*pow(0.7143*(1-0.4*f),1.7);
    if(fabs(aNew-a)/aNew < 0.005){
      a = aNew;
      break;
    }
    a = aNew;
  }
  return a;
}
#include "G4FermiMomentum.hh"
G4double PreCompoundParameters::GetFermiEnergy(const G4Fragment& fragm)
{
  /*  G4double fTmp = fragm.GetA()/(1.07*1.07*1.07*fragm.GetA()*4/3*pi);
      G4FermiMomentum fMom;
      fMom.Init(fragm.GetA(),fragm.GetZ());
      return sqrt(fMom.GetFermiMomentum(fTmp)/MeV*fMom.GetFermiMomentum(fTmp)/MeV + (fTmp=G4NucleiProperties::GetNuclearMass(fragm.GetA(),fragm.GetZ()))*fTmp);*/
  return 35*MeV;
}
