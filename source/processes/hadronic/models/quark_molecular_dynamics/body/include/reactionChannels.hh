#ifndef __REACTIONCHANNELS__
#define __REACTIONCHANNELS__

#include "g4std/vector"
#include "newvector.hh"
#include "genericRead.hh"
#include "EventHandling.hh"

class ParticleBase;
class ParticleType;
class QuantumState;

struct isotop
{
  friend bool operator==(const isotop& a,const isotop& b);
  const ParticleType& pp;
  double charge,mass;
  int anti;
  bool keep;
public:
  isotop(const ParticleType&,double charge_,int anti);
  isotop(const ParticleType&);
  isotop(int C,const ParticleBase&);
};

class CollisionType;

class FunctionType 
{
  double x0,value;
  virtual double f(double x) const = 0;
public:
  double getValue(double x);
};

class ConstantFunction : public FunctionType
{
  double C;
  virtual double f(double) const { return C; }
public:
  ConstantFunction(double c) : C(c) {}
};

class decayMode
{
  friend class CollisionType;
  friend G4std::ostream& operator<<(G4std::ostream& o,decayMode&);

  FunctionType* sigmaPart;
  bool isoSet;
  bool elastic;
  double SumMass;
  selection select;
  vector<isotop*> products;
public:
  decayMode(int C,const vector<ParticleBase*>& P);
  decayMode(G4std::istream&,const vector<ParticleType*>&);
  int N() const { return products.size(); }
  bool isDecomposition() const { return (partialCrossection(1)<0); }
  bool isElastic() const { return elastic; }
  //  bool compareProducts(const vector<ParticleBase*>& P) const;
  double partialCrossection(double s) const { return sigmaPart->getValue(s); }
  void performDecay(const vector<ParticleBase*>& p,double Etot,const Vektor3& beta,const Vektor3& x,bool force = false);
#ifdef  IS_GCC
  decayMode& operator=(const decayMode&) { return *this; }
#endif
};

#endif
