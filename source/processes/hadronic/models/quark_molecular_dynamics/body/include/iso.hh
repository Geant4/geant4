#ifndef __ISO__
#define __ISO__

#include "globals.hh"

class Iso 
{
  static double f(int n,double,double,double,double,double*,double*);
public:
  enum states { ALL=0, LOWEST };
  Iso();
  static void Projections(int,double J,double M,double* j_k,double* m_k,bool* = 0);
  static double chooseMultiplett(double j1,double m1,double j2,double m2,states = ALL);
};

#endif
