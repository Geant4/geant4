#ifndef __METROPOLIS__
#define __METROPOLIS__

#include "g4std/iostream"
#include "newvector.hh"

class Geometry;

class Metropolis
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const Metropolis& M);
protected:
  int N,N_init;
  double T,E_k,F;
  Vektor3* x;
  Vektor3* p;
  int k,l;
  double a,b;
  Geometry& G;
  double E_n(int n) const;
public:
  double kstart;
  Metropolis(Geometry&,int N,double T,double a,double b);
  ~Metropolis();
  void setParticle(const Vektor3& x_,const Vektor3& p_);
  void evaluate(long MAX_ITER);
  double getObservable() const { return F/l; }

  virtual void setX(int n,const Vektor3& x_) {  x[n] = x_; }
  virtual void setP(int n,const Vektor3& p_) {  p[n] = p_; }
  virtual Vektor3 getX(int n) const { return x[n]; }
  virtual Vektor3 getP(int n) const { return p[n]; }
  virtual double Observable() const { return Etot()/N; }
  virtual void doSomething() { }
  virtual double E_part(int i) const = 0;
  virtual double E_int(int i,int j) const = 0;
  virtual double Etot() const = 0;
  virtual void check() {}
};

#endif

