#ifndef __HADRONGAS__
#define __HADRONGAS__

#include "globals.hh"
#include "double.hh"
#include "ParticleType.hh"
#include "bessel.hh"
#include "g4std/iomanip"
#include "g4std/fstream"
#include "newvector.hh"
#include "MathTools.hh"

struct parameter
{
  REAL T,V,pressure;
  REAL mu_q,mu_s;
  Vektor N;
  parameter(REAL T_,REAL muq,REAL mus) 
    : T(T_),V(1.0),mu_q(muq),mu_s(mus),destruct(false),n(0) {}
  parameter(REAL T_,REAL V_,REAL muq,REAL mus) 
    : T(T_),V(V_),mu_q(muq),mu_s(mus),destruct(false),n(0) {}
  parameter(REAL T_,REAL V_,const Vektor& N_) 
    : T(T_),V(V_),N(N_),destruct(false),n(dim(N_)) {}
  ~parameter() {}
  parameter& operator=(const parameter&);
  parameter& operator+=(const parameter&);
  parameter operator*(REAL);

  void init(int);
  void init(int,const Vektor&);
private:
  bool destruct;
  int n;
};

class DistBase
{
public:
  DistBase() {}

  virtual REAL Z(REAL m,REAL T) = 0;
  virtual REAL TdlnZ(REAL m,REAL T) = 0;
  virtual REAL TdTdlnZ(REAL m,REAL T) = 0;
};

class RelBoltzmann : public DistBase
{
public:
  RelBoltzmann() : DistBase() {}

  virtual REAL Z(REAL m,REAL T) { return norm*sqr(m)*T*bessk(2,m/T); }
  virtual REAL TdlnZ(REAL m,REAL T) { REAL mT = m/T; REAL b2 = bessk(2,mT); if ( b2 != 0.0 ) b2 = bessk(1,mT)/b2; else b2 = 1.0; return 3.0+mT*b2; }
  virtual REAL TdTdlnZ(REAL m,REAL T) { REAL mT = m/T; REAL b2 = bessk(2,mT); if ( b2 != 0.0 ) b2 = bessk(1,mT)/b2; else b2 = 1.0; return mT*(mT-4*b2-mT*sqr(b2)); }
public:
  static  REAL norm;
};

class RelBoltzmannStream : public DistBase
{
public:
  RelBoltzmannStream() : DistBase() {}

  virtual REAL Z(REAL m,REAL T) { return 0.5*RelBoltzmann::norm*pow(T,3)*(1+m/T)*exp(-m/T); }
  virtual REAL TdlnZ(REAL m,REAL T) { REAL mT=m/T; return 3.0 + mT*mT/(1+mT); }
  virtual REAL TdTdlnZ(REAL m,REAL T) { REAL mT=m/T; return -mT*(1.0 - 1/sqr(1+mT)); }
};

class gasEnsemble 
{
protected:
public:
  vector<ParticleType*>& List;
  DistBase* K;
public:
  gasEnsemble(DistBase* K_,vector<ParticleType*>& List_) : List(List_),K(K_) {}
  virtual ~gasEnsemble() { delete K; }
  virtual REAL Mu(const ParticleType& h,const parameter&,signed int B) const = 0;
  virtual REAL N(const ParticleType& h,const parameter&,signed int B) const = 0;
  virtual REAL Entropy(const parameter&)  const;
  virtual REAL specificHeat(const parameter&)  const= 0;
  virtual REAL Etot(const parameter&)  const;
  virtual void jacobian(const parameter&,Matrize&) { throw "jacobian:: not defined"; }
  virtual REAL NetN(const ParticleType& h,const parameter& p) const { return N(h,p,1)-N(h,p,-1); }
  virtual REAL SumN(const parameter&,signed int B = 0) const; 
  virtual REAL Pressure(const parameter& p)  const{ return p.T*SumN(p)/p.V; }
  virtual Vektor q_tot(const parameter&) const = 0;
  virtual parameter& initGas(parameter&) { throw "initGas:: not defined"; }
  static const REAL volume;
};

class grandCanonicalEnsemble : public gasEnsemble
{
public:
  grandCanonicalEnsemble( DistBase* K_,vector<ParticleType*>& List_) 
    : gasEnsemble(K_,List_) {}
  REAL specificHeat(const parameter&) const { return 0; }
  //  REAL Entropy(const parameter&) const;
  REAL N(const ParticleType& h,const parameter& p,signed int B = 0) const;
  REAL Mu(const ParticleType& h,const parameter& p,signed int B = 1) const;
  virtual Vektor q_tot(const parameter&) const;
  parameter& initGas(parameter&);
};

class canonicalEnsemble : public gasEnsemble
{
public:
  canonicalEnsemble( DistBase* K_,vector<ParticleType*>& List_) 
    : gasEnsemble(K_,List_) {}
  REAL specificHeat(const parameter&) { return 0; }
  //  REAL Entropy(const parameter&);
  REAL N(const ParticleType& h,const parameter& p,signed int B = 0);
  REAL Mu(const ParticleType& h,const parameter& p,signed int B = 1);
  virtual Vektor q_tot(const parameter&);
  parameter& initGas(parameter& p) { return p; }
};

#endif
