// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                         --- RandPoissonQ ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// M. Fischler    - Implemented new, much faster table-driven algorithm
//		    applicable for mu < 100
//		  - Implemented "quick()" methods, shich are the same as the
//		    new methods for mu < 100 and are a skew-corrected gaussian
//		    approximation for large mu.
// M. Fischler	  - Removed mean=100 from the table-driven set, since it
//		    uses a value just off the end of the table.  (April 2004)
// M Fischler     - put and get to/from streams 12/15/04
// M Fischler     - Utilize RandGaussQ rather than RandGauss, as clearly 
//		    intended by the inclusion of RandGaussQ.h.  Using RandGauss
//		    introduces a subtle trap in that the state of RandPoissonQ
//		    can never be properly captured without also saveing the
//		    state of RandGauss!  RandGaussQ is, on the other hand,
//		    stateless except for the engine used.
// M Fisculer	  - Modified use of wrong engine when shoot (anEngine, mean)
//		    is called.  This flaw was preventing any hope of proper
//		    saving and restoring in the instance cases.
// M Fischler     - fireArray using defaultMean 2/10/05
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// M Fisculer	  - Modified use of shoot (mean) instead of 
//		    shoot(getLocalEngine(), mean) when fire(mean) is called.  
//		    This flaw was causing bad "cross-talk" between modules
//		    in CMS, where one used its own engine, and the other 
//		    used the static generator.  10/18/07
//
// =======================================================================

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Random/Stat.h"
#include "CLHEP/Utility/thread_local.h"
#include <cmath>	// for std::pow()

namespace CLHEP {

std::string RandPoissonQ::name() const {return "RandPoissonQ";}
HepRandomEngine & RandPoissonQ::engine() {return RandPoisson::engine();}

// Initialization of static data:  Note that this is all const static data,
// so that saveEngineStatus properly saves all needed information. 

  // The following MUST MATCH the corresponding values used (in
  // poissonTables.cc) when poissonTables.cdat was created.

const double RandPoissonQ::FIRST_MU = 10;// lowest mu value in table
const double RandPoissonQ::LAST_MU =  95;// highest mu value
const double RandPoissonQ::S = 5;        // Spacing between mu values
const int RandPoissonQ::BELOW = 30;      // Starting point for N is at mu - BELOW
const int RandPoissonQ::ENTRIES = 51;    // Number of entries in each mu row

const double RandPoissonQ::MAXIMUM_POISSON_DEVIATE = 2.0E9;
	// Careful -- this is NOT the maximum number that can be held in 
	// a long.  It actually should be some large number of sigma below
	// that.  

  // Here comes the big (9K bytes) table, kept in a file of 
  // ENTRIES * (FIRST_MU - LAST_MU + 1)/S doubles

static const double poissonTables [ 51 * ( (95-10)/5 + 1 ) ] = {
#include "CLHEP/Random/poissonTables.cdat"
};

//
// Constructors and destructors:
//

RandPoissonQ::~RandPoissonQ() {
}

void RandPoissonQ::setupForDefaultMu() {

  // The following are useful for quick approximation, for large mu
  
  double sig2 = defaultMean * (.9998654 - .08346/defaultMean);
  sigma = std::sqrt(sig2);
	// sigma for the Guassian which approximates the Poisson -- naively
	// sqrt (defaultMean).
	//
	// The multiplier corrects for fact that discretization of the form
	// [gaussian+.5] increases the second moment by a small amount.

  double t = 1./(sig2);

  a2 = t/6 + t*t/324;
  a1 = std::sqrt (1-2*a2*a2*sig2);
  a0 = defaultMean + .5 - sig2 * a2;

  // The formula will be a0 + a1*x + a2*x*x where x has 2nd moment of sigma.
  // The coeffeicients are chosen to match the first THREE moments of the 
  // true Poisson distribution.   
  // 
  // Actually, if the correction for discretization were not needed, then 
  // a2 could be taken one order higher by adding t*t*t/5832.  However,
  // the discretization correction is not perfect, leading to inaccuracy
  // on the order to 1/mu**2, so adding a third term is overkill.  

} // setupForDefaultMu() 


//
// fire, quick, operator(), and shoot methods:
//

long RandPoissonQ::shoot(double xm) {
  return shoot(getTheEngine(), xm);
}

double RandPoissonQ::operator()() {
  return (double) fire();
}

double RandPoissonQ::operator()( double mean ) {
  return (double) fire(mean);
}

long RandPoissonQ::fire(double mean) {
  return shoot(getLocalEngine(), mean);
}

long RandPoissonQ::fire() {
  if ( defaultMean < LAST_MU + S ) {
    return poissonDeviateSmall ( getLocalEngine(), defaultMean );
  } else {
    return poissonDeviateQuick ( getLocalEngine(), a0, a1, a2, sigma );
  }
} // fire()

long RandPoissonQ::shoot(HepRandomEngine* anEngine, double mean) {

  // The following variables, static to this method, apply to the 
  // last time a large mean was supplied; they obviate certain calculations
  // if consecutive calls use the same mean.

  static CLHEP_THREAD_LOCAL double lastLargeMean = -1.;	// Mean from previous shoot 
					// requiring poissonDeviateQuick()
  static CLHEP_THREAD_LOCAL double lastA0;		
  static CLHEP_THREAD_LOCAL double lastA1;		
  static CLHEP_THREAD_LOCAL double lastA2;		
  static CLHEP_THREAD_LOCAL double lastSigma;		

  if ( mean < LAST_MU + S ) {
    return poissonDeviateSmall ( anEngine, mean );
  } else {
    if ( mean != lastLargeMean ) {
      // Compute the coefficients defining the quadratic transformation from a 
      // Gaussian to a Poisson for this mean.  Also save these for next time.
      double sig2 = mean * (.9998654 - .08346/mean);
      lastSigma = std::sqrt(sig2);
      double t = 1./sig2;
      lastA2 = t*(1./6.) + t*t*(1./324.);
      lastA1 = std::sqrt (1-2*lastA2*lastA2*sig2);
      lastA0 = mean + .5 - sig2 * lastA2;
    }
    return poissonDeviateQuick ( anEngine, lastA0, lastA1, lastA2, lastSigma );
  }

} // shoot (anEngine, mean)

void RandPoissonQ::shootArray(const int size, long* vect, double m) {
  for( long* v = vect; v != vect + size; ++v )
    *v = shoot(m);
     // Note: We could test for m > 100, and if it is, precompute a0, a1, a2, 
     // and sigma and call the appropriate form of poissonDeviateQuick.  
     // But since those are cached anyway, not much time would be saved.
}

void RandPoissonQ::fireArray(const int size, long* vect, double m) {
  for( long* v = vect; v != vect + size; ++v )
    *v = fire( m );
}

void RandPoissonQ::fireArray(const int size, long* vect) {
  for( long* v = vect; v != vect + size; ++v )
    *v = fire( defaultMean );
}


// Quick Poisson deviate algorithm used by quick for large mu:

long RandPoissonQ::poissonDeviateQuick ( HepRandomEngine *e, double mu ) {

  // Compute the coefficients defining the quadratic transformation from a 
  // Gaussian to a Poisson:

  double sig2 = mu * (.9998654 - .08346/mu);
  double sig = std::sqrt(sig2);
	// The multiplier corrects for fact that discretization of the form
	// [gaussian+.5] increases the second moment by a small amount.

  double t = 1./sig2;

  double sa2 = t*(1./6.) + t*t*(1./324.);
  double sa1 = std::sqrt (1-2*sa2*sa2*sig2);
  double sa0 = mu + .5 - sig2 * sa2;

  // The formula will be sa0 + sa1*x + sa2*x*x where x has sigma of sq.
  // The coeffeicients are chosen to match the first THREE moments of the 
  // true Poisson distribution.

  return poissonDeviateQuick ( e, sa0, sa1, sa2, sig ); 
} 


long RandPoissonQ::poissonDeviateQuick ( HepRandomEngine *e, 
		double A0, double A1, double A2, double sig) {
//
// Quick Poisson deviate algorithm used by quick for large mu:
//
// The principle:  For very large mu, a poisson distribution can be approximated
// by a gaussian: return the integer part of mu + .5 + g where g is a unit 
// normal.  However, this yelds a miserable approximation at values as
// "large" as 100.  The primary problem is that the poisson distribution is 
// supposed to have a skew of 1/mu**2, and the zero skew of the Guassian 
// leads to errors of order as big as 1/mu**2.
//
// We substitute for the gaussian a quadratic function of that gaussian random.
// The expression looks very nearly like mu + .5 - 1/6 + g + g**2/(6*mu).  
// The small positive quadratic term causes the resulting variate to have 
// a positive skew; the -1/6 constant term is there to correct for this bias 
// in the mean.  By adjusting these two and the linear term, we can match the
// first three moments to high accuracy in 1/mu.
//
// The sigma used is not precisely sqrt(mu) since a rounded-off Gaussian
// has a second moment which is slightly larger than that of the Gaussian.  
// To compensate, sig is multiplied by a factor which is slightly less than 1.

  //  double g = RandGauss::shootQuick( e );   // TEMPORARY MOD:
  double g = RandGaussQ::shoot( e );   // Unit normal
  g *= sig;
  double p = A2*g*g + A1*g + A0;
  if ( p < 0 ) return 0;	// Shouldn't ever possibly happen since 
				// mean should not be less than 100, but
				// we check due to paranoia.
  if ( p > MAXIMUM_POISSON_DEVIATE ) p = MAXIMUM_POISSON_DEVIATE;

  return long(p);

} // poissonDeviateQuick ()



long RandPoissonQ::poissonDeviateSmall (HepRandomEngine * e, double mean) {
  long N1;
  long N2;
  // The following are for later use to form a secondary random s:
  double rRange; 	    // This will hold the interval between cdf for the
			    // computed N1 and cdf for N1+1.
  double rRemainder = 0; // This will hold the length into that interval.

  // Coming in, mean should not be more than LAST_MU + S.  However, we will 
  // be paranoid and test for this:

  if ( mean > LAST_MU + S ) {
    return RandPoisson::shoot(e, mean);
  }

  if (mean <= 0) {
    return 0;			// Perhaps we ought to balk harder here!
  }

  // >>> 1 <<< 
  // 	Generate the first random, which we always will need.

  double r = e->flat();

  // >>> 2 <<< 
  // 	For small mean, below the start of the tables, 
  // 	do the series for cdf directly.  

  // In this case, since we know the series will terminate relatively quickly, 
  // almost alwaye use a precomputed 1/N array without fear of overrunning it.

  static const double oneOverN[50] = 
  {    0,   1.,    1/2.,  1/3.,  1/4.,  1/5.,  1/6.,  1/7.,  1/8.,  1/9., 
   1/10.,  1/11.,  1/12., 1/13., 1/14., 1/15., 1/16., 1/17., 1/18., 1/19., 
   1/20.,  1/21.,  1/22., 1/23., 1/24., 1/25., 1/26., 1/27., 1/28., 1/29.,
   1/30.,  1/31.,  1/32., 1/33., 1/34., 1/35., 1/36., 1/37., 1/38., 1/39., 
   1/40.,  1/41.,  1/42., 1/43., 1/44., 1/45., 1/46., 1/47., 1/48., 1/49. };


  if ( mean < FIRST_MU ) {

    long N = 0;
    double term = std::exp(-mean);
    double cdf = term;

    if ( r < (1 - 1.0E-9) ) {
      //
      // **** This is a normal path: ****
      //
      // Except when r is very close to 1, it is certain that we will exceed r 
      // before the 30-th term in the series, so a simple while loop is OK.
      const double* oneOverNptr = oneOverN;
      while( cdf <= r ) {
        N++ ;  
        oneOverNptr++;
        term *= ( mean * (*oneOverNptr) );
        cdf  += term;
      }
      return N;
      //
      // **** ****
      //
    } else { // r is almost 1...
      // For r very near to 1 we would have to check that we don't fall
      // off the end of the table of 1/N.  Since this is very rare, we just
      // ignore the table and do the identical while loop, using explicit 
      // division.
      double cdf0;
      while ( cdf <= r ) {
        N++ ;  
        term *= ( mean / N );
        cdf0 = cdf;
        cdf  += term;
	if (cdf == cdf0) break; // Can't happen, but just in case...
      }
      return N;
    } // end of if ( r compared to (1 - 1.0E-9) )

  } // End of the code for mean < FIRST_MU

  // >>> 3 <<< 
  // 	Find the row of the tables corresponding to the highest tabulated mu
  //	which is no greater than our actual mean.

  int rowNumber = int((mean - FIRST_MU)/S);
  const double * cdfs = &poissonTables [rowNumber*ENTRIES]; 
  double mu = FIRST_MU + rowNumber*S;
  double deltaMu = mean - mu;
  int Nmin = int(mu - BELOW);
  if (Nmin < 1) Nmin = 1;
  int Nmax = Nmin + (ENTRIES - 1);


  // >>> 4 <<< 
  // 	If r is less that the smallest entry in the row, then 
  //	generate the deviate directly from the series.  

  if ( r < cdfs[0] ) {
  
    // In this case, we are tempted to use the actual mean, and not 
    // generate a second deviate to account for the leftover part mean - mu.
    // That would be an error, generating a distribution with enough excess
    // at Nmin + (mean-mu)/2 to be detectable in 4,000,000 trials.

    // Since this case is very rare (never more than .2% of the r values)
    // and can happen where N will be large (up to 65 for the mu=95 row)
    // we use explicit division so as to avoid having to worry about running
    // out of oneOverN table.

    long N = 0;
    double term = std::exp(-mu);
    double cdf = term;
    double cdf0;

    while(cdf <= r) {
        N++ ;  
        term *= ( mu / N );
        cdf0 = cdf;
        cdf  += term;
	if (cdf == cdf0) break; // Can't happen, but just in case...
    }
    N1 = N;
		// std::cout << r << "   " << N << "   ";
		// DBG_small = true;
    rRange = 0;		// In this case there is always a second r needed

  } // end of small-r case


  // >>> 5 <<< 
  // 	Assuming r lies within the scope of the row for this mu, find the 
  //	largest entry not greater than r.  N1 is the N corresponding to the 
  //	index a.

  else if ( r < cdfs[ENTRIES-1] ) { 		// r is also >= cdfs[0]

    //
    // **** This is the normal code path ****
    //

    int a = 0;                  // Highest value of index such that cdfs[a]
				// is known NOT to be greater than r.
    int b = ENTRIES - 1;	// Lowest value of index such that cdfs[b] is
				// known to exeed r.

    while (b != (a+1) ) {
      int c = (a+b+1)>>1;
      if (r > cdfs[c]) {
        a = c;
      } else {
        b = c;
      }
    }

    N1 = Nmin + a;
    rRange = cdfs[a+1] - cdfs[a];
    rRemainder = r - cdfs[a];

    //
    // **** ****
    //

  } // end of medium-r (normal) case


  // >>> 6 <<< 
  //	If r exceeds the greatest entry in the table for this mu, then start 
  // 	from that cdf, and use the series to compute from there until r is 
  //	exceeded.  

  else { //		  if ( r >= cdfs[ENTRIES-1] ) {
 
    // Here, division must be done explicitly, and we must also protect against
    // roundoff preventing termination.

	// 
	//+++ cdfs[ENTRIES-1] is exp(-mu) sum (mu**m/m! , m=0 to Nmax) 
	//+++ (where Nmax = mu - BELOW + ENTRIES - 1)
	//+++ cdfs[ENTRIES-1]-cdfs[ENTRIES-2] is exp(-mu) mu**(Nmax)/(Nmax)!
	//+++ If the sum up to k-1 <= r < sum up to k, then N = k-1
	//+++ Consider k = Nmax in the above statement:
	//+++ If cdfs[ENTRIES-2] <= r < cdfs[ENTRIES-1], N would be Nmax-1 
	//+++ But here r >= cdfs[ENTRIES-1] so N >= Nmax
	//

	// Erroneous:
	//+++ cdfs[ENTRIES-1] is exp(-mu) sum (mu**m/m! , m=0 to Nmax-1) 
	//+++ cdfs[ENTRIES-1]-cdfs[ENTRIES-2] is exp(-mu) mu**(Nmax-1)/(Nmax-1)!
	//+++ If a sum up to k-1 <= r < sum up to k, then N = k-1
	//+++ So if cdfs[ENTRIES-1] were > r, N would be Nmax-1 (or less)
	//+++ But here r >= cdfs[ENTRIES-1] so N >= Nmax
	//

		//	std::cout << "r = " << r << " mu = " << mu << "\n";
    long N = Nmax -1;
    double cdf = cdfs[ENTRIES-1];
    double term = cdf - cdfs[ENTRIES-2];
    double cdf0;
    while(cdf <= r) {
        N++ ;  
		//	std::cout << "  N " << N << " term " << 
		//	term << " cdf " << cdf << "\n";
        term *= ( mu / N );
	cdf0  = cdf;
        cdf  += term;
	if (cdf == cdf0) break; // If term gets so small cdf stops increasing,
				// terminate using that value of N since we
				// would never reach r.
    }
    N1 = N;
    rRange = 0; 	// We can't validly omit the second true random

	//	    N = Nmax -1;
	//	    cdf = cdfs[ENTRIES-1];
	//	    term = cdf - cdfs[ENTRIES-2];
	//	    for (int isxz=0; isxz < 100; isxz++) {
	//	        N++ ;  
	//	        term *= ( mu / N );
	//		cdf0  = cdf;
	//	        cdf  += term;
	//	    }
	//	    std::cout.precision(20);
	//	    std::cout << "Final sum is " << cdf << "\n";

  } // end of large-r case



  // >>> 7 <<< 
  // 	Form a second random, s, based on the position of r within the range
  //	of this table entry to the next entry.  

  // However, if this range is very small, then we lose too many bits of
  // randomness.  In that situation, we generate a second random for s.

  double s;

  static const double MINRANGE = .01;  	// Sacrifice up to two digits of 
				       	// randomness when using r to produce
					// a second random s.  Leads to up to
					// .09 extra randoms each time.

  if ( rRange > MINRANGE ) {
    //
    // **** This path taken 90% of the time ****
    //
    s = rRemainder / rRange;
  } else {
    s = e->flat();	// extra true random needed about one time in 10.
  }

  // >>> 8 <<< 
  // 	Use the direct summation method to form a second poisson deviate N2 
  //	from deltaMu and s.

  N2 = 0;
  double term = std::exp(-deltaMu);
  double cdf = term;

  if ( s < (1 - 1.0E-10) ) {
      //
      // This is the normal path:
      //
      const double* oneOverNptr = oneOverN;
      while( cdf <= s ) {
        N2++ ;  
        oneOverNptr++;
        term *= ( deltaMu * (*oneOverNptr) );
        cdf  += term;
      }
  } else { // s is almost 1...
      while( cdf <= s ) {
        N2++ ;  
        term *= ( deltaMu / N2 );
        cdf  += term;
      }
  } // end of if ( s compared to (1 - 1.0E-10) )

  // >>> 9 <<< 
  // 	The result is the sum of those two deviates

		// if (DBG_small) {
		//   std::cout << N2 << "   " << N1+N2 << "\n";
		//   DBG_small = false;
		// }

  return N1 + N2;

} // poissonDeviate()

std::ostream & RandPoissonQ::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(a0);
  os << a0 << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(a1);
  os << a1 << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(a2);
  os << a2 << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(sigma);
  os << sigma << " " << t[0] << " " << t[1] << "\n";
  RandPoisson::put(os);
  os.precision(pr);
  return os;
#ifdef REMOVED
  int pr=os.precision(20);
  os << " " << name() << "\n";
  os << a0 << " " << a1 << " " << a2 << "\n";
  os << sigma << "\n";
  RandPoisson::put(os);
  os.precision(pr);
  return os;
#endif
}

std::istream & RandPoissonQ::get ( std::istream & is ) {
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
  if (possibleKeywordInput(is, "Uvec", a0)) {
    std::vector<unsigned long> t(2);
    is >> a0 >> t[0] >> t[1]; a0 = DoubConv::longs2double(t); 
    is >> a1 >> t[0] >> t[1]; a1 = DoubConv::longs2double(t); 
    is >> a2 >> t[0] >> t[1]; a2 = DoubConv::longs2double(t); 
    is >> sigma >> t[0] >> t[1]; sigma = DoubConv::longs2double(t); 
    RandPoisson::get(is);
    return is;
  }
  // is >> a0 encompassed by possibleKeywordInput
  is >> a1 >> a2 >> sigma;
  RandPoisson::get(is);
  return is;
}

}  // namespace CLHEP

