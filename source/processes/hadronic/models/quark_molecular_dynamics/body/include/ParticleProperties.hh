#ifndef __PROPERTIES__
#define __PROPERTIES__

#include "globals.hh"
#include "g4std/iostream"
#include "g4std/vector"
#include "String.hh"
#include "ParticleBase.hh"
#include "ParticleType.hh"

class QuantumProjections
{
protected:
  int c;
  RGB color;
  double iso3,spin3;
public:
  QuantumProjections() : c(1),color(0),iso3(0),spin3(0) {}
  QuantumProjections(int c_,RGB col,double i3,double s3);
  QuantumProjections(const ParticleType&);
  QuantumProjections& operator=(const QuantumProjections& Q);

  double Iso3() const { return iso3; }
  double Spin3() const { return spin3; }
  RGB Color() const { return color; }
  int C() const { return c; }
  void setColor(RGB c_) { color = c_; }
};

class QuantumState : public QuantumProjections,public QuantumNumbers
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const QuantumState& p) { o << p.B() << ": " << p.Spin() << "  " << p.Isospin() << " : " << p.Spin3() << "  " << p.Iso3() << "  "; return o; }
  friend QuantumState operator+(const QuantumState& p1,const QuantumState& p2);
  friend QuantumState addToLowest(const QuantumState& p1,const QuantumState& p2);
  friend QuantumState anti(const QuantumState&);
public:
  QuantumState() {}
  QuantumState(const QuantumNumbers& a,const QuantumProjections& b)
    : QuantumProjections(b),QuantumNumbers(a) {}
  QuantumState(const ParticleBase& p) 
    : QuantumProjections((QuantumProjections&)(ParticleProperties&)p),QuantumNumbers((ParticleType&)p) {}  
  QuantumState& operator=(const QuantumState& x);
};

class ParticleProperties : private QuantumProjections,
			   public virtual ParticleBase
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const ParticleProperties& p) { p.writeOut(o); return o; }

  const ParticleType& Type;
  double mass,lifetime;
	int flag;
	int pdgCode;

protected:
  ParticleProperties(const ParticleType& h,double Emax = 0.0);
  ParticleProperties(const ParticleType& h,const QuantumProjections& q,double Emax = 0.0);
  ParticleProperties(const ParticleProperties& p);
public:
  operator ParticleType() { return Type; }
  void ChargeProjection(double ch);
  //  ParticleProperties anti() const;

  String printCharge() const;
  void setLifetime(double lt = -1.0);

  
  virtual void SetMass(double m) { mass = m; }
  virtual void SetFlag(int f) { flag = f; }
  virtual void SetPDGCode(int pdgc) { pdgCode = pdgc; }
  virtual void SetColor(const RGB& c) { QuantumProjections::setColor(c); }
  virtual void Conjugate() { c = -c; iso3 = -iso3; spin3=-spin3; color=-color; }
  virtual double B() const { return C()*Type.B(); }
  virtual int S() const { return C()*Type.S(); }
  virtual int Charm() const { return C()*Type.Charm(); }
  virtual double Mass() const { return mass; }
  virtual double Isospin() const { return Type.Isospin(); }
  virtual double Spin() const { return Type.Spin(); }

  virtual double Iso3() const { return QuantumProjections::Iso3(); }
  virtual double Spin3() const { return QuantumProjections::Spin3(); }
  virtual RGB Color() const { return QuantumProjections::Color(); }
  virtual int C() const { return QuantumProjections::C(); }

  virtual bool isQuark() const { return Type.isQuark(); }
  virtual bool isDiquark() const { return Type.isDiquark(); }
  virtual bool isGluon() const { return Type.isGluon(); }

  virtual int Flag() const { return flag; }
  virtual String Name() const;

	//
  // PDGCODE
	// depending on how it is set - first guess: constructed from Name()
  //
  virtual int PDGCode() const;
  // 
  // better would be
  //
  //  virtual int PDGCode() const { return pdgCode; }

  virtual double Lifetime() const { return lifetime; }
  virtual QuantumState getProperties() const { return QuantumState(Type,*this); }
  virtual const ParticleType& getType() const { return Type; }
private:
  virtual void writeOut(G4std::ostream& o) const;
  String chargeName() const;
};

#endif
