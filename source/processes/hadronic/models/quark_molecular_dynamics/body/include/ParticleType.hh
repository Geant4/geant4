#ifndef __PARTICLETYPE__
#define __PARTICLETYPE__

#include "ParticleBase.hh"
#include "ReadList.hh"
#include "String.hh"
#include "genericRead.hh"
#include "Tree.hh"
#include "MathTools.hh"
#include "DistributionFunction.hh"
#include "undef.hh"
#include "Pattern.hh"

class G4std::ostream;
class G4std::istream;
class decayMode;
class ParticleDecayModes;
class ParticleProperties;

class RGB
{
  friend RGB operator+(const RGB& c1,const RGB& c2);
  int c;
public:
  enum { WHITE=0,RED,GREEN,BLUE,ANTI_BLUE=-3,ANTI_GREEN,ANTI_RED };
  RGB(int c_ = 0) : c(c_) {}
  RGB operator-() const { return RGB(-c); }
  RGB operator>>(int);
  RGB operator<<(int);
  operator int() const { return c; }
  bool isWhite() const { return c == 0; }
  static RGB pickColor();
};

// 0  -  B
// 1  -  S
// 2  -  spin
// 3  -  isospin
// 4  -  maxcolor
// 5  -  charm
class QuantumNumbers : public Pattern<double>
{
  double peakMass;
public:
  QuantumNumbers() : Pattern<double>(6) {}
  QuantumNumbers(const QuantumNumbers& x) : Pattern<double>(x) {}

  double getIso3() const; 
  double getSpin3() const; 
  RGB getColor(RGB col = 0) const;

  double B() const { return getEntry(0); }
  int S() const { return (int)getEntry(1); }
  int Nq() const { return int(3.0*B()+1e-5)+S(); }
  double Spin() const { return getEntry(2); }
  double Isospin() const { return getEntry(3); }
  int maxColor() const { return (int)getEntry(4); }
  int Charm() const { return (int)getEntry(5); }

  void setB(double B_) { setEntry(0,B_); }
  void setS(int S_) { setEntry(1,S_); }
  void setC(int C_) { setEntry(5,C_); }
  void setSpin(double t3_) { setEntry(2,t3_); }
  void setIsospin(double i3_) { setEntry(3,i3_); }
  void setMaxColor(int c_) { setEntry(4,c_); }
  double PeakMass() const { return peakMass; }
  void setPeakMass(double m_) { peakMass = m_; }

  bool isQuark() const { return fabs(fabs(B())-0.5)<0.5 && fabs(B())<0.5; }
  bool isDiquark() const { return fabs(fabs(B())-0.5)<0.5 && fabs(B())>0.5; }
  bool isGluon() const { return (maxColor() == 8); }

  void writeOut(G4std::ostream& o) const;
};

class ParticleType : public QuantumNumbers,public Knot<ParticleType>
{
  friend class Knot<ParticleType>;
  class MassDist : public GaussIntegration,public SimpleRejection
  {
    const ParticleType& Type;
    double m_max,norm;
    virtual REAL ToBeIntegrated(REAL x);
    virtual REAL f(REAL x) const;
  public:
    MassDist(const ParticleType&,double m);
  };
  friend class MassDist;
  friend class ParticleDecayModes;
  friend G4std::ostream& operator<<(G4std::ostream& o,const ParticleType& x) { x.writeOut(o); return o; }
  friend bool operator==(const ParticleType& x,const ParticleType& y);
  int id,g;
  double minmass,width,BWnorm;
  
  void getProbab(double,vector<double>&,vector<ParticleType*>&,bool = false) const;
  ParticleType();
public:
  ParticleType(G4std::istream&);
  ParticleType& selectType(double mass) const;
  ParticleType& selectType(int,const vector<ParticleBase*>& P,double m) const;
  ParticleType& selectType() const;
  double selectMass(double Mmax) const;
  virtual void ClassInfo(G4std::ostream& o) const { o << " ("<<id << ") "; }
  virtual double getMass(double Emax = 0.0) const;
  double getLifetime() const;
  double getWidth() const { return width; }
  int degeneracy() const { return int(g*(2*Spin()+1)*(2*Isospin()+1)); }
  virtual double isEqualTo(const ParticleType& x) const;

  void writeOut(G4std::ostream& o) const;
  class undefinedParticle {};
private:
};

#endif
