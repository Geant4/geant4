#ifndef __QUARKBOX__
#define __QUARKBOX__

#include "metropolis.hh"
#include "Propagation.hh"
#include "Average.hh"

class QuarkBox : public Metropolis
{
  long cnt;
protected:
  Colour& C;
public:
  QuarkBox(Geometry& G,double T,Colour& C_);
  void setX(int i,const Vektor3& x) { C.List[i]->SetCoordinates(x); }
  void setP(int i,const Vektor3& p) { C.List[i]->SetMomentum(p); }
  Vektor3 getX(int i) const { return C.List[i]->Coordinates(); }
  Vektor3 getP(int i) const { return C.List[i]->Momentum(); }
  double E_part(int i) const { return C.List[i]->E(); }
  double E_int(int i,int j) const { return (i == j) ? 0.0 : C.E_int(i,j); }
  double Etot() const { return C.Etot(); }
  double Eqgp() const;
  double Ehad() const;
  virtual double Observable() const { return Eqgp()/G.getVolume(); }
  double  pressure() const;
  double  phad() const;
  virtual void doSomething();

  Average Eq;
  Average Eh;
  Average Pq;
  Average Ph;
  Average Force;
};

#endif
