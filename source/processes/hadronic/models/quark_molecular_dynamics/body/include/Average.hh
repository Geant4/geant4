#ifndef __AVERAGE__
#define __AVERAGE__

#include "newvector.hh"

class Average
{
  int N;
  double avg;
  double avg2;
public:
  Average() : N(0),avg(0),avg2(0) {}
  void reset() { N=0; avg = 0; avg2 = 0; }
  operator double() { return avg/N; }
  double getError() const { if ( N>1 ) {double x = avg/N; return sqrt((avg2-N*x*x)/(N-1));} }
  void addEntry(double x) { avg += x; avg2 += x*x; ++N; }
  void print(ostream& o) const { o << avg/N << "  " << getError() << endl; }
};

class AverageVector
{
  int N;
  Vektor avg;
  Vektor avg2;
public:
  AverageVector(int n) : N(0),avg(n),avg2(n) {}
  operator Vektor() { return avg/N; }
  Vektor getError() const;
  void addEntry(const Vektor& x);
  void print(ostream& o) const { o << avg/N << "  " << getError() << endl; }
};

#endif
