// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- RandBinomial ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// John Marraffino - Created: 12th May 1998
// M Fischler     - put and get to/from streams 12/10/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
//
// =======================================================================

#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Utility/thread_local.h"
#include <algorithm>	// for min() and max()
#include <cmath>	// for exp()

namespace CLHEP {

std::string RandBinomial::name() const {return "RandBinomial";}
HepRandomEngine & RandBinomial::engine() {return *localEngine;}

RandBinomial::~RandBinomial() {
}

double RandBinomial::shoot( HepRandomEngine *anEngine, long n,
                                                          double p ) {
  return genBinomial( anEngine, n, p );
}

double RandBinomial::shoot( long n, double p ) {
  HepRandomEngine *anEngine = HepRandom::getTheEngine();
  return genBinomial( anEngine, n, p );
}

double RandBinomial::fire( long n, double p ) {
  return genBinomial( localEngine.get(), n, p );
}

void RandBinomial::shootArray( const int size, double* vect,
                            long n, double p )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(n,p);
}

void RandBinomial::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            long n, double p )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(anEngine,n,p);
}

void RandBinomial::fireArray( const int size, double* vect)
{
  for( double* v = vect; v != vect+size; ++v )
    *v = fire(defaultN,defaultP);
}

void RandBinomial::fireArray( const int size, double* vect,
                           long n, double p )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = fire(n,p);
}

/*************************************************************************
 *                                                                       *
 *  StirlingCorrection()                                                 *
 *                                                                       *
 *  Correction term of the Stirling approximation for log(k!)            *
 *  (series in 1/k, or table values for small k)                         *
 *  with long int parameter k                                            *
 *                                                                       *
 *************************************************************************
 *                                                                       *
 * log k! = (k + 1/2)log(k + 1) - (k + 1) + (1/2)log(2Pi) +              *
 *          StirlingCorrection(k + 1)                                    *
 *                                                                       *
 * log k! = (k + 1/2)log(k)     -  k      + (1/2)log(2Pi) +              *
 *          StirlingCorrection(k)                                        *
 *                                                                       *
 *************************************************************************/

static double StirlingCorrection(long int k)
{
  #define   C1               8.33333333333333333e-02     //  +1/12 
  #define   C3              -2.77777777777777778e-03     //  -1/360
  #define   C5               7.93650793650793651e-04     //  +1/1260
  #define   C7              -5.95238095238095238e-04     //  -1/1680

  static const double  c[31] = {   0.0,
			     8.106146679532726e-02, 4.134069595540929e-02,
			     2.767792568499834e-02, 2.079067210376509e-02,
			     1.664469118982119e-02, 1.387612882307075e-02,
			     1.189670994589177e-02, 1.041126526197209e-02,
			     9.255462182712733e-03, 8.330563433362871e-03,
			     7.573675487951841e-03, 6.942840107209530e-03,
			     6.408994188004207e-03, 5.951370112758848e-03,
			     5.554733551962801e-03, 5.207655919609640e-03,
			     4.901395948434738e-03, 4.629153749334029e-03,
			     4.385560249232324e-03, 4.166319691996922e-03,
			     3.967954218640860e-03, 3.787618068444430e-03,
			     3.622960224683090e-03, 3.472021382978770e-03,
			     3.333155636728090e-03, 3.204970228055040e-03,
			     3.086278682608780e-03, 2.976063983550410e-03,
			     2.873449362352470e-03, 2.777674929752690e-03,
  };
  double    r, rr;

  if (k > 30L) {
    r = 1.0 / (double) k;  rr = r * r;
    return(r*(C1 + rr*(C3 + rr*(C5 + rr*C7))));
	}
	else  return(c[k]);
}

double RandBinomial::genBinomial( HepRandomEngine *anEngine, long n, double p )
{
/******************************************************************
 *                                                                *
 *     Binomial-Distribution - Acceptance Rejection/Inversion     *
 *                                                                *
 ******************************************************************
 *                                                                *
 * Acceptance Rejection method combined with Inversion for        *
 * generating Binomial random numbers with parameters             *
 * n (number of trials) and p (probability of success).           *
 * For  min(n*p,n*(1-p)) < 10  the Inversion method is applied:   *
 * The random numbers are generated via sequential search,        *
 * starting at the lowest index k=0. The cumulative probabilities *
 * are avoided by using the technique of chop-down.               *
 * For  min(n*p,n*(1-p)) >= 10  Acceptance Rejection is used:     *
 * The algorithm is based on a hat-function which is uniform in   *
 * the centre region and exponential in the tails.                *
 * A triangular immediate acceptance region in the centre speeds  *
 * up the generation of binomial variates.                        *
 * If candidate k is near the mode, f(k) is computed recursively  *
 * starting at the mode m.                                        *
 * The acceptance test by Stirling's formula is modified          *
 * according to W. Hoermann (1992): The generation of binomial    *
 * random variates, to appear in J. Statist. Comput. Simul.       *
 * If  p < .5  the algorithm is applied to parameters n, p.       *
 * Otherwise p is replaced by 1-p, and k is replaced by n - k.    *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION:    - btpec samples a random number from the binomial *
 *                distribution with parameters n and p  and is    *
 *                valid for  n*min(p,1-p)  >  0.                  *
 * REFERENCE:   - V. Kachitvichyanukul, B.W. Schmeiser (1988):    *
 *                Binomial random variate generation,             *
 *                Communications of the ACM 31, 216-222.          *
 * SUBPROGRAMS: - StirlingCorrection()                            *
 *                            ... Correction term of the Stirling *
 *                                approximation for log(k!)       *
 *                                (series in 1/k or table values  *
 *                                for small k) with long int k    *
 *              - anEngine    ... Pointer to a (0,1)-Uniform      * 
 *                                engine                          *
 *                                                                *
 * Implemented by H. Zechner and P. Busswald, September 1992      *
 ******************************************************************/

#define C1_3     0.33333333333333333
#define C5_8     0.62500000000000000
#define C1_6     0.16666666666666667
#define DMAX_KM  20L

  static CLHEP_THREAD_LOCAL long int      n_last = -1L,  n_prev = -1L;
  static CLHEP_THREAD_LOCAL double        par,np,p0,q,p_last = -1.0, p_prev = -1.0;
  static CLHEP_THREAD_LOCAL long          b,m,nm;
  static CLHEP_THREAD_LOCAL double        pq, rc, ss, xm, xl, xr, ll, lr, c,
				 p1, p2, p3, p4, ch;

  long                 bh,i, K, Km, nK;
  double               f, rm, U, V, X, T, E;

  if (n != n_last || p != p_last)                 // set-up 
	{
	 n_last = n;
	 p_last = p;
	 par=std::min(p,1.0-p);
	 q=1.0-par;
	 np = n*par;

// Check for invalid input values

	 if( np <= 0.0 ) return (-1.0);

	 rm = np + par;
	 m  = (long int) rm;                      // mode, integer 
	 if (np<10)
	{
	 p0=std::exp(n*std::log(q));              // Chop-down
	 bh=(long int)(np+10.0*std::sqrt(np*q));
	 b=std::min(n,bh);
	}
	 else
		 {
	rc = (n + 1.0) * (pq = par / q);          // recurr. relat.
	ss = np * q;                              // variance  
	i  = (long int) (2.195*std::sqrt(ss) - 4.6*q); // i = p1 - 0.5
	xm = m + 0.5;
	xl = (double) (m - i);                    // limit left 
	xr = (double) (m + i + 1L);               // limit right
	f  = (rm - xl) / (rm - xl*par);  ll = f * (1.0 + 0.5*f);
	f  = (xr - rm) / (xr * q);     lr = f * (1.0 + 0.5*f);
	c  = 0.134 + 20.5/(15.3 + (double) m);    // parallelogram
						  // height
	p1 = i + 0.5;
	p2 = p1 * (1.0 + c + c);                  // probabilities
	p3 = p2 + c/ll;                           // of regions 1-4
	p4 = p3 + c/lr;
		 }
  }
  if( np <= 0.0 ) return (-1.0);
  if (np<10)                                      //Inversion Chop-down
	 {
	  double pk;

	  K=0;
	  pk=p0;
	  U=anEngine->flat();
	  while (U>pk)
		{
		 ++K;
		 if (K>b)
			 {
		U=anEngine->flat();
		K=0;
		pk=p0;
			 }
		 else
			 {
		U-=pk;
		pk=(double)(((n-K+1)*par*pk)/(K*q));
			 }
		}
	  return ((p>0.5) ? (double)(n-K):(double)K);
	 }

  for (;;)
	{
	 V = anEngine->flat();
	 if ((U = anEngine->flat() * p4) <= p1)  // triangular region
		{
		 K=(long int) (xm - U + p1*V);
	return ((p>0.5) ? (double)(n-K):(double)K);  // immediate accept
		}
	 if (U <= p2)                                // parallelogram
		{
		 X = xl + (U - p1)/c;
		 if ((V = V*c + 1.0 - std::fabs(xm - X)/p1) >= 1.0)  continue;
		 K = (long int) X;
		}
	 else if (U <= p3)                           // left tail
		{
		 if ((X = xl + std::log(V)/ll) < 0.0)  continue;
		 K = (long int) X;
		 V *= (U - p2) * ll;
		}
	 else                                         // right tail
		{
		 if ((K = (long int) (xr - std::log(V)/lr)) > n)  continue;
		 V *= (U - p3) * lr;
		}

 // acceptance test :  two cases, depending on |K - m|
	 if ((Km = std::labs(K - m)) <= DMAX_KM || Km + Km + 2L >= ss)
	  {

 // computation of p(K) via recurrence relationship from the mode
		f = 1.0;                              // f(m)
		if (m < K)
	 {
	  for (i = m; i < K; )
		{
		if ((f *= (rc / ++i - pq)) < V)  break;  // multiply  f
		}
	 }
		else
	 {
	  for (i = K; i < m; )
		 {
		  if ((V *= (rc / ++i - pq)) > f)  break; // multiply  V
		 }
	 }
		if (V <= f)  break;                       // acceptance test
	 }
  else
	 {

 // lower and upper squeeze tests, based on lower bounds for log p(K)
		V = std::log(V);
		T = - Km * Km / (ss + ss);
		E =  (Km / ss) * ((Km * (Km * C1_3 + C5_8) + C1_6) / ss + 0.5);
		if (V <= T - E)  break;
		if (V <= T + E)
		 {
	if (n != n_prev || par != p_prev)
	 {
	  n_prev = n;
	  p_prev = par;

	  nm = n - m + 1L;
	  ch = xm * std::log((m + 1.0)/(pq * nm)) +
	       StirlingCorrection(m + 1L) + StirlingCorrection(nm);
	 }
	nK = n - K + 1L;

 // computation of log f(K) via Stirling's formula
 // final acceptance-rejection test
	if (V <= ch + (n + 1.0)*std::log((double) nm / (double) nK) +
                 (K + 0.5)*std::log(nK * pq / (K + 1.0)) -
                 StirlingCorrection(K + 1L) - StirlingCorrection(nK))  break;
		}
	 }
  }
  return ((p>0.5) ? (double)(n-K):(double)K);
}

std::ostream & RandBinomial::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultP);
  os << defaultN << " " << defaultP << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
}

std::istream & RandBinomial::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  if (possibleKeywordInput(is, "Uvec", defaultN)) {
    std::vector<unsigned long> t(2);
    is >> defaultN >> defaultP;
    is >> t[0] >> t[1]; defaultP = DoubConv::longs2double(t); 
    return is;
  }
  // is >> defaultN encompassed by possibleKeywordInput
  is >> defaultP;
  return is;
}


}  // namespace CLHEP
