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
  P=sqrt(E*E-1.0) ;
  U=Z/137.0;
  S=sqrt(1.0-U*U) - 1.;
  Y = 2*PI*U*E/P;
  A1 = U*U*E*E + P*P/4.;
  A2 = fabs(Y/(1-exp(-Y)));
  F2 = pow(A1,S) * A2; 
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
  S=sqrt(1.0-U*U) - 1.;
  G4double F1 = 1E-10;
  for (G4int i = 1; i<=100 ; i++) {
    E = G4double(i)*EE + 1.;
    P=sqrt(E*E-1.0) ;
    Y = 2*PI*U*E/P;
    A1 = U*U*E*E + P*P/4.;
    A2 = fabs(Y/(1-exp(-Y)));
    F2 = pow(A1,S) * A2; 
    if (F2 > F1) F1 = F2;
  }		   
  return F1;
}



















