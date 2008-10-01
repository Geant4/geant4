// usage demonstrated in unit/GNUmakefile
#include "unit/G4Unit.hh"
#include <iostream>
#include <math.h>

using namespace std;

// Based on INUCL pin.f and output:
//  IN,IT,KZ,S,PM,AMPI =    3   1   1  0.392E+00  0.146E+00  0.140E+00
//  L,PHASE =    2 -0.794E+01
//  L,PHASE =    5 -0.165E+01
//  L,PHASE =    6  0.823E+01
//       CALL PHFUN(3,1,0,S31,PM,S,AMPI) ! S 31 PHASE
//       CALL PHFUN(3,1,1,P31,PM,S,AMPI) ! P 31 PHASE
//       CALL PHFUN(3,3,1,P33,PM,S,AMPI) ! P 33 PHASE
// we have data for unit testing:
// phfun(3, 1, 0, 0.146, 0.392, 0.140)= -7.94
// phfun(3, 1, 1, 0.146, 0.392, 0.140)= -1.65
// phfun(3, 3, 1, 0.146, 0.392, 0.140)=  8.23

//       SUBROUTINE PHFUN(I,J,K,PHASE,PM,S,AMPI)
// *    CALC. TANG(DELTA) FOR DIFF. T,L
//       PARAMETER (PI=3.1415927, TWOPI=6.2831854)
//       DIMENSION X(8),S0(8),XK0(8),G0(8),B(8),C(8),D(8)
//       DATA X/0.44,0.31,0.61,0.23,0.22,0.99,0.54,0.43/
//       DATA S0/1.55,1.655,1.435,1.815,1.85,1.233,1.525,1.67/
//       DATA XK0/0.477,0.55,0.393,0.656,0.678,0.228,0.459,0.56/
//       DATA G0/0.105,0.17,0.23,0.255,0.20,0.116,0.125,0.155/
//       DATA B/16.8,-11.2,-5.71,-1.31,-2.91,11.4,0.109,0.112/
//       DATA C/-35.4,-30.7,25.4,1.22,3.45,-15.4,-0.031,-0.27/
//       DATA D/27.,21.,-29.,-0.4,-1.5,7.2,0.003,0.19/
// *     L DEFINITION
//       IF(I.EQ.1.AND.J.EQ.1.AND.K.EQ.0) L=1  ! S11
//       IF(I.EQ.3.AND.J.EQ.1.AND.K.EQ.0) L=2  ! S31
//       IF(I.EQ.1.AND.J.EQ.1.AND.K.EQ.1) L=3  ! P11
//       IF(I.EQ.1.AND.J.EQ.3.AND.K.EQ.1) L=4  ! P13
//       IF(I.EQ.3.AND.J.EQ.1.AND.K.EQ.1) L=5 ! P31
//       IF(I.EQ.3.AND.J.EQ.3.AND.K.EQ.1) L=6 ! P33
//       IF(I.EQ.1.AND.J.EQ.3.AND.K.EQ.2) L=7  ! D13
//       IF(I.EQ.1.AND.J.EQ.5.AND.K.EQ.2) L=8  ! D15
//       XR=X(L)
//       S0R=S0(L)
//       XK0R=XK0(L)
//       G0R=G0(L)
//       BR=0.01*B(L)
//       CR=0.001*C(L)
//       DR=0.0001*D(L)
// *   O'K, NOW ONLY PHASE CALC. REMAINES
//       XX=PM/AMPI
// C      PRINT '('' XX= '',E10.3)',XX
//       TAND=XX**(2*K+1)*(BR+CR*XX**2+DR*XX**4)
//      *          +
//      *     XR*((PM/XK0R)**(2*K+1))*G0R*S0R/(S0R**2-S)
// C      PRINT '('' TAND = '',E10.3)',TAND
//       PHASE=ATAN(TAND)
//       IF(L.EQ.6.AND.PHASE.LE.0.) PHASE=PI+PHASE
//       PHASE=2.*PHASE
// C!!!!!!!!!!!!!!!!!!!!!!!!  -->C
//       PRINT '('' L,PHASE = '',I4,(1X,E10.3))',L,PHASE*90./PI
//       END

double phfun(int i, int j, int k, double pm, double s, double ampi) {
// Calculate tang(delta) for diff. t,l 
// translated from INCUL PHFUN by aatos.heikkinen@cern.ch

	const double   x[8] = {0.44, 0.31, 0.61, 0.23, 0.22, 0.99, 0.54, 0.43};

	const double  s0[8] = {1.550, 1.655, 1.435, 1.815, 1.850, 1.233, 1.525, 1.670};
	const double xk0[8] = {0.477, 0.550, 0.393, 0.656, 0.678, 0.228, 0.459, 0.560};
	const double  g0[8] = {0.105, 0.170, 0.230, 0.255, 0.200, 0.116, 0.125, 0.155};

	const double   b[8] = {  16.8, -11.2,  -5.71, -1.31, -2.91,  11.4,  0.109,  0.112};
	const double   c[8] = { -35.4, -30.7,  25.40,  1.22,  3.45, -15.4, -0.031, -0.270};
	const double   d[8] = {  27.0,  21.0, -29.00, -0.40, -1.50,   7.2,  0.003,  0.190};

// l definition
	int l = 0;
	if(i==1 && j==1 && k==0) l=1;  // S11
	if(i==3 && j==1 && k==0) l=2;  // S31
	if(i==1 && j==1 && k==1) l=3;  // P11
	if(i==1 && j==3 && k==1) l=4;  // P13
	if(i==3 && j==1 && k==1) l=5;  // P31
	if(i==3 && j==3 && k==1) l=6;  // P33
	if(i==1 && j==3 && k==2) l=7;  // D13
	if(i==1 && j==5 && k==2) l=8;  // D15

	double   xr =   x[l];
	double  s0r =  s0[l];
	double xk0r = xk0[l];
	double  g0r =  g0[l];

	double br = 0.01   * b[l];
	double cr = 0.001  * c[l];
	double dr = 0.0001 * d[l];

// Phase calculation remaines
	double xx = pm/ampi;

//       TAND=XX**(2*K+1)*(BR+CR*XX**2+DR*XX**4)
//      *          +
//      *     XR*((PM/XK0R)**(2*K+1))*G0R*S0R/(S0R**2-S)

	double tand = pow(xx, 2*k+1) * (br + cr * pow(xx, 2) + dr * pow(xx,4)) + 
		xr * pow(pm/xk0r, 2*k+1) * g0r * s0r / (pow(s0r, 2) - s);

	double phase = atan(tand);

	double pi=3.1415927; 
	if(l==6 && phase <= 0.0) phase = pi + phase;
	phase = 2.0 * phase;

	std::cout << "l = "<< l << ", xx = "<< xx << " tand = "  << tand 
		  << " phase = " << phase*90/pi << std::endl;

	return phase;
}


testSuite(IllustrationSuite)
namespace {

	testCase(FloatingPoint, IllustrationSuite)
	{
		using namespace BERT;
		using namespace std;

		float expectedValue = 1.234f;
		float errorLimit = 0.001f;
		float value = 1.03f + 0.204f;
		// test if two float values were close enough to each other
		assertTrue(closeValue(expectedValue, value, errorLimit));

		float phaseS31 = phfun(3, 1, 0, 0.146, 0.392, 0.140);  
 
		cout << endl << "phase" << endl; 
		cout << "S31 " << phaseS31 << endl;

		assertTrue(closeValue(-7.94, phaseS31, errorLimit)); // phfun(3, 1, 0, 0.146, 0.392, 0.140)= -7.94

		float phaseP31 = phfun(3, 1, 1, 0.146, 0.392, 0.140);
		float phaseP33 = phfun(3, 3, 1, 0.146, 0.392, 0.140);

		assertTrue(closeValue(expectedValue, value, 0.5));
//      assertTrue(closeValue(1.0, 2.0, 0.1)); // AH::: fails    

		// specific information when an assertion failed: BERT::closeValueInfo()
		// NOTE: closeValueInfo() uses operator << to convert parameters to string description
		// be sure a proper operator << is provided for each parameter
		// it is already defined for float, double and long double
		assertTrue(closeValueInfo(expectedValue, value, errorLimit));
		assertTrue(closeValueInfo(1.66666666, 5.0 / 3.0, 0.001));

		// there is also a function for equality test: BERT::equalValueInfo()
		// NOTE: proper operator << for every parameter is also needed
		const int expectedInt = 8;
		assertTrue(equalValueInfo(expectedInt, 0x01 << 3));
	}

}// namespace
