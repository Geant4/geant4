#ifndef __NEWBINNING__
#define __NEWBINNING__

#include "globals.hh"
#include "newvector.hh"

class Binning
{
  friend ostream& operator<<(ostream& o,const Binning& B) { B._print(o); return o; }
  int n;
protected:
  long Ntot;
  bool locked;
  double xmin,xmax,dx;
  Vektor value;
  Vektor err;
  int getIndex(double x) const;
  void _print(ostream&) const;
  virtual double normalize(double) const { return Ntot*dx; }
public:
  Binning() : n(0) {}
  Binning(int n_,double xmin_,double xmax_);
  Binning(const Binning&);
  virtual ~Binning() {}

  Vektor getValues() const { return value; }
  Vektor getRange() const;

  void init(int n_,double xmin_,double xmax_);
  double getMax() const;
  void addEntry(double x,double y = 1.0);
  double getEntry(int i) const;
  void lock();
  Binning& operator*=(double x);
  Binning& operator/=(double x);
  Binning operator/(const Binning& B) const;
  Binning operator*(double x) const;
  Binning operator+(const Binning& B) const;
  Binning operator-(const Binning& B) const;
  virtual double x(int i) const { return xmin+(i-0.5)*dx; }

  class outOfBounds {};
  class notLocked {};
};

#endif
