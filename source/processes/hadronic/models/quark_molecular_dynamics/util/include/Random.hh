#ifndef __RANDOM__
#define __RANDOM__

#include <stdlib.h>

class G4std::ostream;

//<etymology>
// rand_gen is a contraction for <b>rand</b>om number <b>gen</b>erator.
//</etymology>
//<motivation>
// Provide an easy-to-use Random Generator for floating point as well as integer type variables.
//</motivation>
//<synopsis>
// We used a 48-bit random number generated provided in the C-library of many UNIX systems. If your machine does not support this kind of random number generator please define the functions <src>drand48</src> and <src>srand48</src> used in the code with the respective functions of your random generator. Further, the system time is used in order to give a "random" seed. Therefore the function <src>gettimeofday</src> and the <src>struct timeval</src> is used that is provided in <src>sys/time.h</src> on AIX machines. Please ask your system administrator for an adequat replacement on your system and include this definition in file <src>function.C</src>.
//</synopsis>
class rand_gen
{
// Print seeding of the random generator. This is important for event reproduction.
  friend G4std::ostream& operator<<(G4std::ostream&,rand_gen&);
public:
// Initializes random generator with seed given by argument. 0 initializes the
// random generator by using the actual system time.
  rand_gen(long int = 0);
  void Set(long seed);
// Gives random number of type REAL (floating point) in the range 0 &lt= x &lt 1.
  double operator()() const { return drand48(); }
// Returns a random number between 0 and <src>max</src>.
  double operator()(double max) const { return max*drand48(); }
// Returns a random number in the range <src>min</src>&lt= x &lt<src>max</src>.
  double operator()(double min,double max) const { return min+(max-min)*drand48(); }
  int operator()(int min,int max) const { return int(min+(max+1-min)*drand48()); }
  int randInt(double x) { int n = (int)x; return ( drand48() < x - n) ? n+1 : n; }
  static rand_gen Random;
private:
  long int rand_seed;
};


#endif

