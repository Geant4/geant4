#ifndef __PARTICLEBASE__
#define __PARTICLEBASE__

#include "newvector.hh"
#include "globals.hh"
#include "g4std/iostream"

class String;
class QuantumState;
class QuantumProjections;
class ParticleType;
class RGB;

class ParticleBase
{
  friend G4std::ostream& operator<<(G4std::ostream&,const ParticleBase&);
public:
  virtual ~ParticleBase() {}
  virtual int C() const = 0;
  virtual double B() const = 0;
  virtual int S() const = 0;
  virtual int Charm() const = 0;
  virtual RGB Color() const = 0;
  virtual double Mass() const = 0;
  virtual int Flag() const = 0;
  virtual double Isospin() const = 0;
  virtual double Spin() const = 0;
  virtual double Iso3() const = 0;
  virtual double Spin3() const = 0;
  virtual String Name() const = 0;
  virtual double Time() const = 0;
  virtual double Lifetime() const = 0;
  virtual QuantumState getProperties() const = 0;
  virtual const ParticleType& getType() const = 0;

  virtual void SetMass(double m) = 0;
  virtual void SetFlag(int flag) = 0;
  virtual void setLifetime(double) = 0;
  virtual void SetColor(const RGB& c) = 0;
  virtual void Conjugate() = 0;
  virtual Vektor3 Coordinates() const = 0;
  virtual double Coordinates(int k) const = 0;
  virtual double& Coordinates(int k) = 0;
  virtual void SetCoordinates(const Vektor3& x) = 0;
  virtual void SetCoordinates4(const Vektor4& x) = 0;
  virtual Vektor3 Momentum() const = 0;
  virtual double Momentum(int k) const = 0;
  virtual double& Momentum(int k) = 0;
  virtual void SetMomentum(const Vektor3& x) = 0;
  virtual void SetMomentum4(const Vektor4& x) = 0;
  virtual ParticleBase* makeClone() const = 0;

  double TransverseMass() const; 
  double TransverseMomentum() const; 

  virtual bool isQuark() const = 0;
  virtual bool isDiquark() const = 0;
  virtual bool isGluon() const = 0;

  virtual Vektor3 Force() const { return Vektor3(0,0,0); }
  void Lorentz(int j,double beta);
  void Lorentz(const Vektor3& beta);
  double Charge() const { return Iso3()+0.5*(B()+S()+Charm()); }
  double E() const  { return sqrt(Momentum()*Momentum()+Mass()*Mass()); }
protected:
  virtual void announceEvent() const {} 
  virtual bool noPotentials() const = 0;
  virtual double V(double r) const = 0;
  virtual double dVdr(double r) const = 0;
  virtual Vektor3 Vprime(const ParticleBase&) = 0;
};

//ParticleBase* makeParticle(const String&);
ParticleBase* makeParticle(const ParticleType&);
ParticleBase* makeParticle(const ParticleType&,const QuantumProjections&,double mass = -1);
ParticleBase* makeParticle(const ParticleType&,const QuantumProjections&,const Vektor3& p,const Vektor3& x,double mass = -1);

#endif
