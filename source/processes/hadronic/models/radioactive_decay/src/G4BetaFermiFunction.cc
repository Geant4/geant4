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
#include "globals.hh"
#include "G4BetaFermiFunction.hh"

const G4double G4BetaFermiFunction::PI=3.14159;

//////////////////////////////////////////////////////////////////
//
// calculate the Fermi Function foe energy E0
//
G4double G4BetaFermiFunction::GetFF( const G4double E0)
{
  G4double A1, A2;
  G4double P, U, S, Y;
  G4double F2;
  G4double E = E0+1.;  
  P=std::sqrt(E*E-1.0) ;
  U=Z/137.0;
  S=std::sqrt(1.0-U*U) - 1.;
  Y = 2*PI*U*E/P;
  A1 = U*U*E*E + P*P/4.;
  A2 = std::fabs(Y/(1-std::exp(-Y)));
  F2 = std::pow(A1,S) * A2; 
  return F2;
}

//////////////////////////////////////////////////////////////////
//
//  calculate the Fermi normalization factor 
//  here E0 is the end point energy of the beta decay
//
G4double G4BetaFermiFunction::GetFFN(const G4double E0)
{

  G4double A1, A2;
  G4double P, U, S, Y;
  G4double F2,E;
  G4double EE = E0/100.;
  U=Z/137.0;
  S=std::sqrt(1.0-U*U) - 1.;
  G4double F1 = 1E-10;
  for (G4int i = 1; i<=100 ; i++) {
    E = G4double(i)*EE + 1.;
    P=std::sqrt(E*E-1.0) ;
    Y = 2*PI*U*E/P;
    A1 = U*U*E*E + P*P/4.;
    A2 = std::fabs(Y/(1-std::exp(-Y)));
    F2 = std::pow(A1,S) * A2; 
    if (F2 > F1) F1 = F2;
  }		   
  return F1;
}



















