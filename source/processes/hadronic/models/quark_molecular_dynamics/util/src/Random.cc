#include "Random.hh"
#include <iostream.h>

rand_gen rand_gen::Random;

#ifdef __PC__
   #include "mytime.hh"
#define throw exit(0); 
#else
   #include <sys/time.h>
#endif

void rand_gen::Set(long seed) 
{
  if (!seed)
    {
      struct timeval tv;
      gettimeofday(&tv,0);
      rand_seed = tv.tv_sec + tv.tv_usec;
      srand48((unsigned int)(rand_seed));
    }
  else {
    srand48((unsigned int)seed);
    rand_seed = seed;
  }
}

rand_gen::rand_gen(long int seed) : rand_seed(seed)
{
  Set(seed);
}

ostream& operator<<(ostream& o,rand_gen& r)
{
	o << "Random Seed = " << r.rand_seed << endl;
	return o;
}












