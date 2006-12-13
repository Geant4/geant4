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
#include <math.h>
#include "G4FermiFunction.hh"

//////////////////////////////////////////////////////////////////
//
//
double G4FermiFunction::GetFF( const double E0)
{
  double A1, A2;
  double P, U, S, Y;
  double F2;
  double E = E0+1.;  
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
//
double G4FermiFunction::GetFFN(const double E0)
{

  double A1, A2;
  double P, U, S, Y;
  double F2,E;
  double EE = E0/100.;
  U=Z/137.0;
  S=std::sqrt(1.0-U*U) - 1.;
  double F1 = 1E-10;
  for (int i = 1; i<=100 ; i++) {
    E = double(i)*EE + 1.;
    P=std::sqrt(E*E-1.0) ;
    Y = 2*PI*U*E/P;
    A1 = U*U*E*E + P*P/4.;
    A2 = std::fabs(Y/(1-std::exp(-Y)));
    F2 = std::pow(A1,S) * A2; 
    if (F2 > F1) F1 = F2;
  }		   
  return F1;
}


//////////////////////////////////////////////////////////////////
//
//
double G4FermiFunction::Gamma(const double X) 
{
  static double   PI = 4.0*std::atan(1.0);

  //              COEFFICIENTS FOR MINIMAX 
  //              APPROXIMATION TO GAMMA(X),
  //              2.0 .LE. X .LE. 3.0
  static double P[5] = {-51.49952, 80.05398,-201.4659,-1.889439, 9.895546};
  static double Q[4] = {130.5263, -303.5898, 26.84174, -19.52375};
  //                                  APPROXIMATION TO LN(GAMMA(X)),
  //                               12.0 .LE. X
  static double P4[3] = {.9189385, .8333332E-01, -.2770927E-02};
  //
  static int IEND = 4 , IEND1 = 3,  IEND2 = 2;
  static double  XINF = 1.7E+38;
  //  GAMMA(XMIN) .APPROX. XINF
  //      GAMMA(BIG1) .APPROX. XINF
  static double             XMIN = 5.8775E-39;
  static double             BIG1 = 34.844;


  //               SPECIFICATIONS FOR LOCAL VARIABLES
  double      XSIGN,Y,T,R,A,TOP,DEN,B;
  int         I, J;

  bool        MFLAG;

  double GAMMA;
  //                                  FIRST EXECUTABLE STATEMENT

  MFLAG = false;
  T = X;
  if (std::fabs(T) <= XMIN) 
    {
      GAMMA = XINF;
      if (T <= 0.0) GAMMA = -XINF;
      //cout << GAMMA<< endl;
      return GAMMA;
    }
  else if (std::fabs(T) >= BIG1) 
    {
      GAMMA = XINF;
      //      cout << GAMMA <<endl;
      return GAMMA;
    }
  else if (T <= 0.0) 
    {
      //                                  ARGUMENT IS NEGATIVE
      MFLAG = true;
      T = -T;
      R = (T);
      XSIGN = 1.0;
      if ( std::fmod(R,2.0) == 0.0) XSIGN = -1.;
      R = T-R;
      if(R == 0.0) {
	GAMMA = XINF;
	if (XSIGN == -1.0) GAMMA = -XINF;
	//cout << GAMMA << endl;
	return GAMMA;
      } 
      else
	{
	  //                                  ARGUMENT IS NOT A NEGATIVE INTEGER
	  R = PI/ std::sin(R*PI) * XSIGN;
	  T += 1.0 ;}
    } 
  else if (T <= 12.0) 
    {
      //                              EVALUATE APPROXIMATION FOR GAMMA(T)
      //                              T .GT. XMIN
      I = int(T);
      A = 1.0;
      if (I == 0) {
	//                                  0.0 .LT. T .LT. 1.0
	A = A/(T*(T+1.0));
	T = T+2.0;}
      else if (I ==1) 
	{      
	  //                                  1.0 .LE. T .LT. 2.0
	  A /= T;   
	  T += 1.0;}
      else if (I > 2) 
	{      
	  //                                  3.0 .LE. T .LE. 12.0
	  for ( J= 3; J<= I; J++ )
	    {
	      T -= 1.0;
	      A *= T;
	    }
	}
      //                                  2.0 .LE. T .LE. 3.0
      TOP = P[IEND1]*T+P[IEND];
      DEN = T+Q[IEND1];
      for ( J=0; J<=IEND2; J++) 
	{   
	  TOP = TOP*T+P[J];
	  DEN = DEN*T+Q[J];}
      Y = (TOP/DEN)*A;
      if (MFLAG) Y = R/Y;
      GAMMA = Y ;
      //      cout << X <<"  "<< GAMMA<<endl;
      return GAMMA;
    }
  else 
    {
      //                                  T .GT. 12.0
      TOP = std::log(T);
      TOP = T*(TOP-1.0)-.5*TOP;
      T = 1.0/T ;
      B = T*T ;
      Y = (P4[2]*B+P4[1])*T+P4[0]+TOP;
      Y = std::exp(Y);
      if (MFLAG) Y = R/Y;
      GAMMA = Y;
      //      cout << X<< "  "<< GAMMA<<endl;
      return GAMMA;
    }
}


















