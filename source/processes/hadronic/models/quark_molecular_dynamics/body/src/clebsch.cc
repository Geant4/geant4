#include <math.h>
#include "clebsch.hh"

ClebschGordan ClebschGordan::Clebsch;

ClebschGordan::ClebschGordan()
{
  logfak[0] = 0.;
  for (int i = 1; i <= 100; ++i) {
    logfak[i] = logfak[i - 1] + log((double)i);
  }
} 

double ClebschGordan::operator()(double j1, double m1, double j2, double m2, double j3,double m3)
{
    double sig = 1.0;
    if ( int(-j1+j2-m3) % 2 )
      sig = -1;
    return sig*sqrt(2.0*j3+1.0)*w3j(j1,j2,j3,m1,m2,-m3);
} 

double ClebschGordan::w3j(double j1, double j2, double j3, double m1, double m2, double m3)
{
/*  This program calculates the 3-j wigner symbols according to the */
/*  representation of A. Lindner. */

/*  Reference: */
/*  A. Lindner, Drehimpulse in der Quantenmechanik, Teubner 1984, P.39 */

/*  Written by Ch. Hofmann, 16.7.92, last changes: 22.4.93 */

/* ====================================================================== 
*/
/*  Program input */
/*  Program returns */
/*  Global variables */
/*  Program variables */
/*  Start of calculation */
/*  Evaluation due to equivalence with Regge symbol */

    double sigma = j1 + j2 + j3;
    double n[9];
    n[0] = -(j1) + j2 + j3;
    n[3] = j1 - j2 + j3;
    n[6] = j1 + j2 - j3;
    n[1] = j1 - m1;
    n[4] = j2 - m2;
    n[7] = j3 - m3;
    n[2] = j1 + m1;
    n[5] = j2 + m2;
    n[8] = j3 + m3;
    for (int i = 1; i <= 3; ++i) {
      for (int j = 1; j <= 3; ++j) {
	if (int(n[i + j * 3 - 4]) < 0) {
	  return 0.0;
	}
      }
      double sum1 = n[i - 1] + n[i + 2] + n[i + 5];
      double sum2 = n[i * 3 - 3] + n[i * 3 - 2] + n[i * 3 - 1];
      if (int(sum1) != int(sigma) || int(sum2)!= int(sigma)) 
	return 0.0;
    }
    int imin = 1;
    int jmin = 1;
    int signum = 1;
    double minimal = n[0];
/*  Looking for the smallest N( i, j ) */
    {for (int i = 1; i <= 3; ++i) {
	for (int j = 1; j <= 3; ++j) {
	    if (n[i + j * 3 - 4] < minimal) {
		minimal = n[i + j * 3 - 4];
		imin = i;
		jmin = j;
	    }
	}
    }}
    signum = 1;
    if (imin > 1) {
	for (int j = 1; j <= 3; ++j) {
	    double dummy = n[j * 3 - 3];
	    n[j * 3 - 3] = n[imin + j * 3 - 4];
	    n[imin + j * 3 - 4] = dummy;
	}
	if ( int(sigma) % 2 )
	  signum = -1;
    }
    if (jmin > 1) {
	for (int i = 1; i <= 3; ++i) {
	    double dummy = n[i - 1];
	    n[i - 1] = n[i + jmin * 3 - 4];
	    n[i + jmin * 3 - 4] = dummy;
	}
	if ( int(sigma) % 2 )
	  signum = -signum;
    }
    double r1 = n[0];
    double r2 = n[3];
    double r3 = n[6];
    double r4 = n[1];
    double r5 = n[4];
    double r6 = n[7];
    double r7 = n[2];
    double r8 = n[5];
    double r9 = n[8];
    double lf_r1__ = logfak[int(r1)];
    double lf_r2__ = logfak[int(r2)];
    double lf_r3__ = logfak[int(r3)];
    double lf_r4__ = logfak[int(r4)];
    double lf_r5__ = logfak[int(r5)];
    double lf_r6__ = logfak[int(r6)];
    double lf_r7__ = logfak[int(r7)];
    double lf_r8__ = logfak[int(r8)];
    double lf_r9__ = logfak[int(r9)];
    double lf_sigma__ = logfak[int(sigma)+1];
    double hlp1 = (lf_r2__ + lf_r3__ + lf_r4__ + lf_r7__ - lf_sigma__ - lf_r1__ -   
	    lf_r5__ - lf_r9__ - lf_r6__ - lf_r8__) / 2.;
    double pre = exp(hlp1) * ((int(r6 + r8) % 2) ? -1.0 : 1.0);
    double hlp2 = lf_r6__ - logfak[int(r6 - r1)] + lf_r8__ - 
	    logfak[int(r8 - r1)];
    double s[101];
    s[0] = exp(hlp2);
    double summe = s[0];
    for (int in = 1; in <= int(r1); ++in) {
	double dn = (double) in;
	s[in] = -s[in - 1] * (r1 + 1. - dn) * (r5 + 1. - dn) * (r9 + 1. - dn) 
		/ dn / (r6 - r1 + dn) / (r8 - r1 + dn);
	summe += s[in];
    }
    double ret_val = pre * summe * signum;
    return ret_val;
}


