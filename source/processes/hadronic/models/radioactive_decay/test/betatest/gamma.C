#include <stdlib.h>
#include <iostream.h>

#include <math.h>

int main ()

//double Gamma(G4double const X) 
{
  double X;
  
  //               SPECIFICATIONS FOR LOCAL VARIABLES

  static double   PI = 4.0*atan(1.0);

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

  double      XSIGN,Y,T,R,A,TOP,DEN,B;
  int         I, J;

  int        MFLAG;

  double GAMMA;
  //                                  FIRST EXECUTABLE STATEMENT

  for (int k = 0; k<=100; k++){
 
    X = k/10.;

  MFLAG = 0;
  T = X;
  if (fabs(T) <= XMIN) 
    {
      GAMMA = XINF;
      if (T <= 0.0) GAMMA = -XINF;
      cout << GAMMA<< endl;
    }
  else if (fabs(T) >= BIG1) 
    {
      GAMMA = XINF;
      cout << GAMMA <<endl;
      //      return GAMMA;
    }
  else if (T <= 0.0) 
    {
      //                                  ARGUMENT IS NEGATIVE
      MFLAG = 1;
      T = -T;
      R = (T);
      XSIGN = 1.0;
      if ( fmod(R,2.0) == 0.0) XSIGN = -1.;
      R = T-R;
      if(R == 0.0) {
	GAMMA = XINF;
	if (XSIGN == -1.0) GAMMA = -XINF;
	cout << GAMMA << endl;
	//return GAMMA;
      } 
      else
	{
	  //                                  ARGUMENT IS NOT A NEGATIVE INTEGER
	  R = PI/ sin(R*PI) * XSIGN;
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
      cout << X <<"  "<< GAMMA<<endl;
    }
  else 
    {
      //                                  T .GT. 12.0
      TOP = log(T);
      TOP = T*(TOP-1.0)-.5*TOP;
      T = 1.0/T ;
      B = T*T ;
      Y = (P4[2]*B+P4[1])*T+P4[0]+TOP;
      Y = exp(Y);
      if (MFLAG) Y = R/Y;
      GAMMA = Y;
      cout << X<< "  "<< GAMMA<<endl;
    }
  } 

}
