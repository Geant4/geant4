#ifndef __MATH_NEU__
#define __MATH_NEU__
#ifdef IS_GCC
#pragma interface
#endif

#include <math.h>

template<class t>
t sqr(t x) { return x*x; }

template<class t>
t sign(t x) { return (x>0 ? 1 : ((x<0) ? -1 : 0)); }

#endif
