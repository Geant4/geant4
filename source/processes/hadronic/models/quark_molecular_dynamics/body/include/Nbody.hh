#ifndef __NBODY__
#define __NBODY__

#include "g4std/vector"
#include <string.h>
#include "Random.hh"
#include "newvector.hh"
#include "g4std/iostream"
#include <math.h>
#include "String.hh"
#include "ParticleBase.hh"

class Particle;

class Nbody
{
  friend class Particle;

  double time;
  double h;
  Vektor array;
  Vektor forces;
  void handleCollisions();
protected:
  void add(Particle*,int = -1);
  void sub(Particle*);
  void sub(int i);
  virtual Vektor& function(Vektor& f); 
  void checkDecay();

  bool firstCall;
  //  int N,Npart;
  virtual inline double Etot();
  virtual inline Vektor3 dHdp(int i);
  virtual inline Vektor3 dHdx(int i) { return Vektor3(); }
public:
  int N,Npart;
  vector<Particle*> List;
  Nbody(double = 1e-6);
  virtual ~Nbody();
  virtual void checkRange() {}
  inline virtual Vektor3 dr(const Vektor3& x1,const Vektor3& x2) const;
  inline virtual double distance(int i,int j) const;
  virtual void one_step();
  virtual void print(G4std::ostream& o);

  void setTime(double t) { time = t; }
  double Time() const { return time; }
  double TimeStep() const { return h; }
};

class Particle : public virtual ParticleBase
{
  friend class Nbody;
  int offset;
  double time;
  double rho;
protected:
  double L,norm;
  Particle();
  Particle(int pos);
  Particle(const Vektor3& mom,const Vektor3& x,int pos);
  virtual Vektor3 Vprime(const ParticleBase&);
  double gauss(const Vektor3& x) const;
public:
  virtual ~Particle();
  //  Particle& operator=(const ParticleProperties& p) { ParticleProperties::set(p); return *this; }
  virtual double getParam(int) const { return 0.0; }
  virtual void setParam(int,double) { }
  virtual void reset() {}
  virtual void refresh() {}
  virtual double Time() const { return time; }
  virtual inline Vektor3 Coordinates() const;
  virtual inline double Coordinates(int k) const;
  virtual inline double& Coordinates(int k);
  virtual void SetCoordinates(const Vektor3& x);
  virtual void SetCoordinates(const Vektor3& x,double);
  virtual void SetCoordinates4(const Vektor4& x);
  virtual inline Vektor3 Momentum() const;
  virtual inline double Momentum(int k) const;
  virtual inline double& Momentum(int k);
  virtual void SetMomentum(const Vektor3& x);
  virtual void SetMomentum4(const Vektor4& x);
  virtual inline Vektor3 Force() const;
  inline Vektor3 Velocity() const;

  static Nbody* Soup;
  static double rho_0;
};

#include "Nbody.icc"

#endif
