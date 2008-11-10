// usage demonstrated in unit/GNUmakefile
#include "unit/G4Unit.hh"
#include <iostream>
#include <math.h>

typedef double G4double;
typedef int G4int;

G4int G4VerbosityLevel=0;  // global verbosity 
using namespace std;
G4double pi    = 3.1415927; 
G4double twopi = 6.2831854;
G4double inDeg = pi/90;


G4double phfun(const G4int i, const G4int j, const G4int k, 
               const G4double pm, const G4double s, const G4double ampi) {
        // Calculate tan(delta) for different pi kinetic energes and l.
        // Data and fitting function from R. Rowe et al. Phys. Review C 18(1) pp.584-589 
        // translated from inuclcnew.f function PHFUN 
        // ampi pion mass

// *     PARAMETERS: P(3) - SCM MOMENTA, IN - PION TYPE, IT - NUCL., *
// *     KZ - RECH. KEY, S=ESCM(PI)**2 (GEV), PM - !P SCM!,          *
// *     AMPI - PION MASS
                                            
	const G4double   x[8] = {0.44, 0.31, 0.61, 0.23, 0.22, 0.99, 0.54, 0.43};

	const G4double  s0[8] = {1.550, 1.655, 1.435, 1.815, 1.850, 1.233, 1.525, 1.670};
	const G4double xk0[8] = {0.477, 0.550, 0.393, 0.656, 0.678, 0.228, 0.459, 0.560};
	const G4double  g0[8] = {0.105, 0.170, 0.230, 0.255, 0.200, 0.116, 0.125, 0.155};

	const G4double   b[8] = {  16.8, -11.2,  -5.71, -1.31, -2.91,  11.4,  0.109,  0.112};
	const G4double   c[8] = { -35.4, -30.7,  25.40,  1.22,  3.45, -15.4, -0.031, -0.270};
	const G4double   d[8] = {  27.0,  21.0, -29.00, -0.40, -1.50,   7.2,  0.003,  0.190};


	G4int l = 0; // l definition
	if(i==1 && j==1 && k==0) l=1;  // S11
	if(i==3 && j==1 && k==0) l=2;  // S31
	if(i==1 && j==1 && k==1) l=3;  // P11
	if(i==1 && j==3 && k==1) l=4;  // P13
	if(i==3 && j==1 && k==1) l=5;  // P31
	if(i==3 && j==3 && k==1) l=6;  // P33
	if(i==1 && j==3 && k==2) l=7;  // D13
	if(i==1 && j==5 && k==2) l=8;  // D15

	G4int lc = l-1; // cpp indexing 
	G4double   xr =   x[lc];
	G4double  s0r =  s0[lc];
	G4double xk0r = xk0[lc];
	G4double  g0r =  g0[lc];

	G4double br = 0.01   * b[lc];
	G4double cr = 0.001  * c[lc];
	G4double dr = 0.0001 * d[lc];

        // Phase calculation remains
	G4double xx = pm / ampi;

	G4double tand = pow(xx, 2*k+1) * 
		(br + cr * pow(xx, 2) + dr * pow(xx, 4)) + 
		xr * pow(pm/xk0r, 2*k+1) * g0r * s0r / (pow(s0r, 2) - s);
	G4double phase = atan(tand);

	if(l==6 && phase <= 0.0) phase = pi + phase;
	phase = 2.0 * phase;
        if (G4VerbosityLevel>0) {
		cout << "\t\t\t l = "      << l << ", xx = "<< xx << "  pi mass = " << ampi
		     << ", phase = " << phase / inDeg << " deg" << endl;
	}
	return phase;
}

//       SUBROUTINE TETNEW(P,IN,IT,KZ,S,PM,AMPI,IRET)
// *******************************************************************
// *     CALCULATE GOOD TETA DISTRIBUTIONS FOR PI - N INTERACTIONS   *
// *     BELOW 400. MEV ENERGY                                       *
// *     ACCORDING TO ROWE,SOLOMON,LANDAW                            *
// *     N.STEPANOV, 13.4.93                                         *
// *                                                                 *
// *     PARAMETERS: P(3) - SCM MOMENTA, IN - PION TYPE, IT - NUCL., *
// *     KZ - RECH. KEY, S=ESCM(PI)**2 (GEV), PM - !P SCM!,          *
// *     AMPI - PION MASS                                            *
// *******************************************************************
//       PARAMETER (PI=3.1415927, TWOPI=6.2831854)
//       DIMENSION P(3)
//       IRET=0    ! IF -10 RETURN TO NO "DELTA"
// C!!!!!!!!!!!!!!!!!!!!!!!!   -->CC
// C     PRINT '('' IN,IT,KZ,S,PM,AMPI = '',3I4,3(1X,E10.3))',
// C    *           IN,IT,KZ,S,PM,AMPI
// *     D(SIGMA)/D(TETA) CALC.FOR T = 3/2
//       CALL PHFUN(3,1,0,S31,PM,S,AMPI) ! S 31 PHASE
//       CALL PHFUN(3,1,1,P31,PM,S,AMPI) ! P 31 PHASE
//       CALL PHFUN(3,3,1,P33,PM,S,AMPI) ! P 33 PHASE
// *     RE (A3/2)
//       X32=COS(S31)-1.
//       Y32=2.*(COS(P33)-1.)+COS(P31)-1.  ! COS TERM
// *     IM (A3/2)
//       Z32=SIN(S31)
//       W32=2.*SIN(P33)+SIN(P31)      ! COS TERM
// *     RE (B3/2)
//       REB32=COS(P33)-COS(P31)
// *     IM (B3/2)
//       B32IM=SIN(P33)-SIN(P31)
// *
// *     NOW CALC. 1/2 TERM , IF NEEDED AND ALL OTHER COEFF.
//       IF((IN.EQ.3.AND.IT.EQ.1).OR.(IN.EQ.5.AND.IT.EQ.2)) THEN
// *       PURELY 3/2 CASE
//         D32=1.
//         D12=0.
//        ELSE ! REAL CALC. OF 1/2 TERM NEEDED
//         CALL PHFUN(1,1,0,S11,PM,S,AMPI) ! S 11 PHASE
//         CALL PHFUN(1,1,1,P11,PM,S,AMPI) ! P 11 PHASE
//         CALL PHFUN(1,3,1,P13,PM,S,AMPI) ! P 13 PHASE
// *       2.5.93 D13,D15 ARE RESERVED FOR BETTER TIMES
// C!!!!!!!!!!!!!!!!!!!!!!!!!  -->CC
// C       CALL PHFUN(1,3,2,D13,PM,S,AMPI) ! D 13 PHASE
// C       CALL PHFUN(1,5,2,D15,PM,S,AMPI) ! D 15 PHASE
// *     RE (A1/2)
//       X12=COS(S11)-1.
//       Y12=2.*(COS(P13)-1.)+COS(P11)-1.  ! COS TERM
// *     IM (A1/2)
//       Z12=SIN(S11)
//       W12=2.*SIN(P13)+SIN(P11)      ! COS TERM
// *     RE (B1/2)
//       REB12=COS(P13)-COS(P11)
// *     IM (B3/2)
//       B12IM=SIN(P13)-SIN(P11)
// *
// *      NOW RIGHT COEFF., WEIGHTING 3/2 AND 1/2 TERMS
//         IF(KZ.EQ.1) THEN ! ELASTIC SCAT.
//           IF((IN.EQ.5.AND.IT.EQ.1).OR.(IN.EQ.3.AND.IT.EQ.2)) THEN
// *       PI- P  OR PI+ N  CASE
//             D32=1./SQRT(3.)
//             D12=SQRT(2./3.)
//            ELSE ! PI 0 CASE
//             D32=SQRT(2./3.)
//             D12=SQRT(1./3.)
//            ENDIF
//          ELSE ! RECHARGERING
//           D32=1.  !    SQRT(2.)/3
//           D12=-1.
//          ENDIF
//        ENDIF
// C!!!!!!!!!!!!!!!!!!!!!!!!!!    -->C
// C     PRINT '('' D32,D12 = '',2(E11.4,2X))',D32,D12
// *    O'K NOW REAL WEIGHTED COEFF.
//       XX=D32*X32+D12*X12
//       YY=D32*Y32+D12*Y12
//       ZZ=D32*Z32+D12*Z12
//       WW=D32*W32+D12*W12
//       BB=(D32*REB32+D12*REB12)**2+(D32*B32IM+D12*B12IM)**2
//       A=YY**2+WW**2-BB  ! COS**2 TERM
//       B=2.*(XX*YY+ZZ*WW)  ! COS TERM
//       C=XX**2+ZZ**2+BB  ! CONST. TERM
// ***************************
// C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   -->CC
// C     PRINT '('' A,B,C,E,F = '',5(E10.3,1X))',A,B,C,E,F
// C     STOP                !  FOR A,B,C CHECK
// ***************************
// *  O'K NOW TETA , ACCORDING TO A*COS**2+B*COS+C DISTRIBUTION
// *      MAX VALUE OF DENSITY DEFINITION
// *     TEST OF DISTRIBUTION
//       DET=B**2-4.*A*C
//       IF(DET.GT.0.) THEN
//         X1=0.5*(SQRT(DET)-B)/A
//         X2=0.5*(-SQRT(DET)-B)/A
//         IF(ABS(X1).LT.1..OR.ABS(X2).LT.1.) THEN ! BAD DISTRIBUTION
//           IRET=-10
//           PRINT '('' TETNEW: SMTH. WRONGE '')'
//           PRINT '('' IN,IT,KZ,S,PM,AMPI = '',3I4,3(1X,E10.3))',
//      *             IN,IT,KZ,S,PM,AMPI
//           PRINT '('' A,B,C = '',3(E10.3,1X))',A,B,C
//           RETURN
//          ENDIF
//        ENDIF
//       IF(A.GE.0.) THEN ! MAX. ALWAYS AT THE -1. OR +1.
//         F1M=A-B+C
//         F2M=A+B+C
//         FMAX=MAX(F2M,F1M)
//        ELSE ! HERE ARE 3 POSSIBLE POINT OF MAX.
//         XGM=-0.5*B/A  ! POINT OF GLOBAL MAX.
//         IF(XGM.GT.1.) THEN !  LOC. MAX. AT 1.
//           FMAX=A+B+C
//          ELSEIF(XGM.LT.-1.) THEN ! LOC. MAX. AT -1.
//           FMAX=A-B+C
//          ELSE ! GLOBAL ONE
//           FMAX=C-0.25*B**2/A
//          ENDIF
//        ENDIF
// C     PRINT '('' FMAX = '',E10.3)',FMAX
//       IF(FMAX.LE.0.) THEN
//         PRINT '('' TETNEW: SMTH. WRONGE : FMAX = '',E10.3)',FMAX
//         PRINT '('' IN,IT,KZ,S,PM,AMPI = '',3I4,3(1X,E10.3))',
//      *             IN,IT,KZ,S,PM,AMPI
//         IRET=-10.
//         RETURN
//        ENDIF
// C...........................................................RANDOM
// 1     S=RNDM(-1)*FMAX
// C...........................................................RANDOM
//       COST=1.-2.*RNDM(-2)  ! COS OF C.M.S ANGLE
//       F=A*COST**2+B*COST+C
//       IF(F.LT.S) GOTO 1
// *     O'K NOW CALC. P(1-3)
// C     CALL HF1(100,COST,1.)
// C...........................................................RANDOM
//       FI=TWOPI*RNDM(-3)
//       PT=PM*SQRT(1.-COST**2)
//       P(1)=PT*COS(FI)
//       P(2)=PT*SIN(FI)
//       P(3)=PM*COST
//       END


//       SUBROUTINE TETNEW(P,IN,IT,KZ,S,PM,AMPI,IRET)
G4int tetnew(G4int in, G4int it, G4int kz, G4double s, G4double pm, G4double ampi) {
// Calculate good teta distributions for pi-N interactions   
// below 400. MeV energy. According to R. Rowe et al. Phys. Review C 18(1) pp.584-589 
// translated from INCUL TETNEW

// Parameters: p(3) - scm momenta, in - pion type, it - nucl., 
// kz - rech. key, s=pow(escm(pi), 2) (GeV), pm - !p scm!,          
// ampi - pion mass                                            
	G4double p[3];
	G4int iret=0; // if -10 return to no "delta"
        if (G4VerbosityLevel>0) {
		cout << "\t\t\t in = "<< in << ", it = "<< it << " kz = "  << kz 
		     << " s = " << s  << " pm = " << pm << " ampi = " << ampi << endl;
	}
// d(sigma) / d(teta) calculation for t = 3/2
	G4double s31 = phfun(3, 1, 0, pm, s, ampi); // s31 phase
	G4double p31 = phfun(3, 1, 1, pm, s, ampi); // p31 phase
	G4double p33 = phfun(3, 3, 1, pm, s, ampi); // p33 phase

	G4double x32 = cos(s31)-1.0; // Re (a3/2)
	G4double y32 = 2.0 * (cos(p33) - 1.0) + cos(p31) - 1.0;  // cos-term

	G4double z32 = sin(s31); // Im (a3/2)
	G4double w32 = 2.0*sin(p33) + sin(p31);      // cos-term

	G4double reb32 = cos(p33) - cos(p31); // Re (b3/2)

	G4double b32im = sin(p33) - sin(p31); // Im (b3/2)

        G4double d32;
        G4double d12;
        G4double x12;
	G4double y12;
	G4double z12;
	G4double w12;
        G4double reb12;
	G4double b12im;
	G4double fmax;
        if (G4VerbosityLevel>0) {
		cout << "\t\t\t tetnew: s31 = "  << s31*90/pi << " deg." << 
			" p31 = " << p31 * 90/pi << " deg." << 
			" p33 = " << p33 * 90/pi << " deg." << endl;
	}
        // now calculate 1/2 term , if needed and all other coefficients
	if((in==3 && it==1) || (in==5 && it==2)) { // purely 3/2 case
		d32 = 1.0;
		d12 = 0.0;
	} else 	{ // real calculation of 1/2 term needed

		G4double s11 = phfun(1,1,0,pm,s,ampi); // s11 phase
		G4double p11 = phfun(1,1,1,pm,s,ampi); // p11 phase
		G4double p13 = phfun(1,3,1,pm,s,ampi); // p13 phase
//G4double d13 = phfun(1,3,2,pm,s,ampi); // d13 phase (d13 and d15 are reserved for future use)
//G4double d15 = phfun(1,5,2,pm,s,ampi); // d15 phase
		// Re (a1/2)
		x12 = cos(s11) - 1.0;
		y12 = 2.0 * (cos(p13) - 1.0) + cos(p11) - 1.0; // cos -term

		// Im (a1/2)
		z12 = sin(s11);
		w12 = 2.0 * sin(p13) + sin(p11); // cos -term

                // Re (b1/2)
		reb12 = cos(p13) - cos(p11);

                // Im (b3/2)
		b12im = sin(p13) - sin(p11);

// Now right coefficients, weighting 3/2 and 1/2 terms

		if (kz == 1) { // elastic scattering
			if ((in == 5 && it == 1) || (in == 3 && it == 2)) { // pi- p or pi+ n case

				G4double d32=1.0 / sqrt(3.0);
				G4double d12 = sqrt(2.0 / 3.0);
			} else { // pi0 case
				G4double d32 = sqrt(2.0 / 3.0);
				G4double d12 = sqrt(1.0 / 3.0);
			}
		} else { // rechargering
			G4double d32 = 1; // sqrt(2)/3;
			G4double d12 = -1.0;

		}
	}
	if (G4VerbosityLevel>0) {
		cout << "d32, d12 = " << d32 << ", " << d12 << endl;
	}

// o'k now real weighted coefficient
	G4double xx = d32 * x32 + d12 * x12;
	G4double yy = d32 * y32 + d12 * y12;
	G4double zz = d32 * z32 + d12 * z12;
	G4double ww = d32 * w32 + d12 * w12;
	G4double bb = pow(d32 * reb32 + d12 * reb12, 2) + pow(d32 * b32im + d12 * b12im, 2);

	G4double a = pow(yy, 2) + pow(ww, 2) - bb; // pow(cos, 2) term
	G4double b = 2.0 * (xx * yy + zz * ww);    //        cos  term
	G4double c = pow(xx, 2) + pow(zz, 2) + bb; //   constant  term
	if (G4VerbosityLevel>0) {
		cout << "a, b, c = " << a << ", " << ", "<< b << ", " << c << endl;
	}

// now theta sampling according to a * pow{cos, 2) + b * cos + c distribution
// max value of density definition
// test of distribution
	G4double det = pow(b, 2) - 4.0 * a * c;

	if (det > 0.0) {
		G4double x1 = 0.5 * (sqrt(det)  - b) / a;
		G4double x2 = 0.5 * (-sqrt(det) - b) / a;

		if (fabs(x1) < 1.0 || fabs(x2) < 1.0) { // bad distribution 
			G4int iret = -10;

//			if (G4VerbosityLevel>0) {
			cout << "ERROR: Bertini cascade / tetnew: something wronge" << endl;
			cout <<  "in, it, kz, s, pm, ampi = " << in << ", " << it << ", " << 
				kz << ", " << s << ", " << pm << ", " << ampi << endl;
			cout << "a, b, c = " << a << ", " << ", "<< b << ", " << c << endl;
//			}
			return iret;
		}
	}
	if(a >= 0.) { // maximum always at the -1.0 or +1.0
		G4double f1m = a - b + c;
		G4double f2m = a + b + c;
		fmax = max(f2m, f1m);

	} else { // here are 3 possible point of maximum
		G4double xgm = -0.5 * b / a;  // point of global maximum
		if(xgm > 1.0) { //  local maximum at 1.0
			fmax = a + b + c;
		} else if(xgm < -1.0) { // local maximum at -1.0
			fmax = a - b + c;
		} else { // global one
			G4double fmax = c - 0.25 * pow(b, 2) / a;
		}
	}


	if (G4VerbosityLevel>0) {
		cout << "fmax = " << fmax << endl;
	}
	if (fmax <= 0.0) {
// if (G4VerbosityLevel>0) {
		cout << "ERROR Bertini cascade / tetnew" << endl;
		cout << "something wrong : fmax " << fmax << endl;
		cout << "in, it, kz, s, pm, ampi = " << in << ", " << it << 
			", " << kz << ", " << s << ", " << pm << ", " << ampi << endl;
//	}
		iret = -10;
		return iret;
	}
//AH:::: HOT SPOT >>>
// C...........................................................RANDOM
// 1     S=RNDM(-1)*FMAX
// C...........................................................RANDOM
//       COST=1.-2.*RNDM(-2)  ! COS OF C.M.S ANGLE
//       F=A*COST**2+B*COST+C
//       IF(F.LT.S) GOTO 1
//      O'K NOW CALC. P(1-3)
// C     CALL HF1(100,COST,1.)
// C...........................................................RANDOM
//       FI=TWOPI*RNDM(-3)
//       PT=PM*SQRT(1.-COST**2)
//       P(1)=PT*COS(FI)
//       P(2)=PT*SIN(FI)
//       P(3)=PM*COST
	return iret;  //AH::: return P[]
}



// C... TWO USE CASES EXIST IN INUCLCNEW.F
// C...1ST
//         IF(IN.LT.3.OR.EKL.GT.0.4) THEN ! TETA CASE
//           CALL TETA(P17,IS,PP,EKL,KZ)
//          ELSE
//           CALL TETNEW(P17,IN,IT,KZ,SSQ,PP,AM(IN),IRET)
//           IF(IRET.LT.0) CALL TETA(P17,IS,PP,EKL,KZ)
//          ENDIF
// C...2ND
//       IF(IN0.LT.3.OR.EK.GT.400.) THEN ! TETA CASE
//         CALL TETA(P17,IN0,PSCM,EK,KZ)
//        ELSE
//         CALL TETNEW(P17,IN0,1,KZ,SSQ,PSCM,AMIN,IRET)
//         IF(IRET.LT.0) CALL TETA(P17,IN0,PSCM,EK,KZ)
//        ENDIF
// C...END USE CASES

//  IN,IT,KZ,S,PM,AMPI =    3   2   2  0.107E+00  0.172E+00  0.140E+00
//  L,PHASE =    2 -0.101E+02
//  L,PHASE =    5 -0.255E+01
//  L,PHASE =    6  0.121E+02
//  L,PHASE =    1  0.906E+01
//  L,PHASE =    3 -0.220E+01
//  L,PHASE =    4 -0.117E+01

testSuite(BertiniSuite) // Bertini unit testing
namespace {
	testCase(Test, BertiniSuite) {
		using namespace BERT;
		using namespace std;

		G4VerbosityLevel = 0;

		if (G4VerbosityLevel>0) cout << "\n testCase / test:\n";

		G4double tolerance = 0.2;
		if (G4VerbosityLevel>0) {
			cout << "\t test\n";
		}
	};

	testCase(PionPhase_phfun, BertiniSuite) {
		using namespace BERT;
		using namespace std;
		G4VerbosityLevel=1;

		if (G4VerbosityLevel>0) cout << "\n testCase / PionPhase_phfun:\n";
 
                G4double piMinusMass = 0.13957; // Bertini masses  
                G4double piZeroMass  = 0.13498; 
		G4double errMax      = 0.2; // tolerance deg.

		float phaseS31 = phfun(3, 1, 0, 0.146, 0.392, 0.140);  
		if (G4VerbosityLevel>0) {
			cout << endl << "\t phase S31 = " << phaseS31 / inDeg << " deg." << endl;
		}

		assertTrue(closeValue(-7.94, phfun(3, 1, 0, 0.146, 0.392, 0.140) / inDeg, errMax)); 
		if (G4VerbosityLevel>0) {
			cout << "\t phase S31: sensitivity analysis:"<< endl;
			cout << " \t \t1. variation of pi mass (tolerance errMax = " << 
				errMax << " deg\n)" << endl;
		}



        	assertTrue(
			closeValue(-7.94, phfun(3, 1, 0, 0.146, 0.392, piMinusMass) / inDeg, errMax)); 

	
// Must fail
		cout << "\t WRONG MASS SHOULD FAIL:"<< endl;
        	assertTrue(
			closeValue(-7.94, phfun(3, 1, 0, 0.146, 0.392, piZeroMass ) / inDeg, errMax));     

//		assertTrue(closeValue(-1.65, (phfun(3, 1, 1, 0.146, 0.392, 0.140) / inDeg), errMax)); 
//		assertTrue(closeValue( 8.23, (phfun(3, 3, 1, 0.146, 0.392, 0.140) / inDeg), errMax)); 
	};


	testCase(PionPhase_tetnew, BertiniSuite)   {
		using namespace BERT;
		using namespace std;

		G4VerbosityLevel=0;

		if (G4VerbosityLevel>0) cout << "\n testCase / PionPhase_tetnew:\n";

		if (G4VerbosityLevel>0) {
			cout << "\t tetnew = " << tetnew(3, 2, 2, 0.107, 0.172, 0.140) << endl;
		}
	
	}
}// namespace
