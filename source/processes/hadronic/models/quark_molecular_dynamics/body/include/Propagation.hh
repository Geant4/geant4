#ifndef __PROPAGATION__
#define __PROPAGATION__

#include <vector.h>
#include "newvector.hh"
#include "Nbody.hh"
#include "ParticleType.hh"
#include "ParticleProperties.hh"
#include "reactionChannels.hh"
#include "PotentialBase.hh"
#include "MathTools.hh"
#include "Geometry.hh"
#include "double.hh"

struct connect
{
  friend ostream& operator<<(ostream& o,connect& x) { o << x.pointer << "  " << x.dist<< "  " << x.max << "  "<< x.force << "  "; return o;}
  int pointer;
  double dist,max,force;
  double t1,t2;
public:
  connect() : dist(1e30),pointer(0),max(0),force(0),t1(0),t2(0) {}
  void set(int i,double d) { pointer = i; dist = d;}
  void set(int i,double d,double m) { pointer = i; dist = d; max = m;}
  void reset() { dist=1e30; }
};

class Quark;

class Colour : public Nbody,protected InverseFunctionWithNewton
{
protected:
  enum { MAX_TRIES = 20 };
  friend class Quark;
  static double colors[7][7];
  static double colors_close[7][7];
  static double* FlavorFission;
  static double Schwinger(double kap);
  ParticleType& TunnelRate(double kap);
 
  Vektor3 f,r_k;
  double lf,sigma;
  int kk;
  REAL ToBeInverted(REAL x) const;
  REAL Derivative(REAL x) const;
  vector<int> eraseList;
public:
  enum criterion { ALL, NO_DIQUARK };
  virtual void init(double) {}
  Colour(double h = 0.01);
  ~Colour();
  inline double Etot();
  inline double E(int i,int j);
  inline double E(int i);
  inline double E_int(int i,int j);
  inline double V(double r,double s = 1) const { return Pot->V(r,s); }
  inline double V_inv(double r,double s = 1) const { return Pot->V_inv(r,s); }
  inline double V_prime(double r,double s = 1) const { return Pot->V_prime(r,s); }
  inline Vektor3 dHdx(int i);
  Vektor3 Ptot();
  void one_step();
  virtual void print(ostream& o);
  void Correlation();
  void clusters(int i);
  void decomposition(int i);
  void Fission();
  void setup();
  void decompose(ParticleBase*);
  void decomposeAll();
  double field(const Vektor3& r);
  void writeField(ostream&,double,double,double,double,double,double);
  static ParticleType& selectQuark(criterion = ALL);
  static double kappa;
  static double alpha;
  static double alpha_s;
  static double A;
  static double B,a,R;
  static double shift;
  static double shift2;
  static double FissionRate;
  static bool decay;
  static bool allowClustering;
  static bool allowDecomposition;
  static bool allowDecay;
  static bool directHadrons;
  static bool removeWhenUndefined;
  static double deconfined;  // Units: [1]  (force relative to kappa)
  static double minDist;  // Units: [1]  (corresponds to factors of av. hadron radius)
  static double decompDist;  // Units: [fm]
  static PotentialBase* Pot;
  static vector<ParticleType*> Quarks;
  static void setQuarks(vector<ParticleType*>&);
  int Nquark;
  double Factor(double r,int i,int j,int s = 1) const;
private:
  static double colorProfile(double r) { return 1.0; }
  double a0,r0,E0;
};

class Radiation : public Colour
{
  Geometry& G;
public:
  Radiation (Geometry& g,double h = 0.01) 
    : Colour(h),G(g) {}
  virtual void init(double);
  virtual void checkRange();
  Vektor3 dHdx(int i) { return List[i]->Color() ? Vektor3(0,0,kappa) : Vektor3(0,0,0); }
  inline double Etot();
  inline double V(double) const;
  inline double V_inv(double) const;
  inline double V_prime(double) const;
  void one_step();
  void FindCorrelations();
};

class inTheBox : public Colour
{
  Geometry& G;
public:
  inTheBox(Geometry& g,double h = 0.01) 
    : Colour(h),G(g) {}
  virtual void checkRange();
  virtual Vektor3 dr(const Vektor3&,const Vektor3&) const;
};

class RealParticle : protected Particle,protected ParticleProperties,
		     public virtual ParticleBase
{
protected:
  RealParticle(const ParticleType& h,int pos)
    : Particle(pos),ParticleProperties(h) { setLifetime(); }
  RealParticle(const ParticleType& h,const QuantumProjections& qp,int pos)
    : Particle(pos),ParticleProperties(h,qp) { setLifetime(); }
  RealParticle(const ParticleType& h,const Vektor3& p,const Vektor3& x,int pos)
    : Particle(p,x,pos),ParticleProperties(h) { setLifetime(); }
  RealParticle(const ParticleType& h,const QuantumProjections& qp,const Vektor3& p,const Vektor3& x,int pos)
    : Particle(p,x,pos),ParticleProperties(h,qp) { setLifetime(); }
};

class Quark : public RealParticle
{
  double path,traj_dt;
  double FissionProbab;
  double justFiss;
public:
  connect next[3];
  void reset();
  Quark(const ParticleType& h) : RealParticle(h,0) { reset(); L*=pow(3.0,0.33333); norm /= 3.0; }
  Quark(const ParticleType& h,const QuantumProjections& qp) 
    : RealParticle(h,qp,0) { reset(); L*=pow(3.0,0.33333); norm /= 3.0; }
  Quark(const ParticleType& h,const Vektor3& p,const Vektor3& x)
    : RealParticle(h,p,x,0) { reset();L*=pow(3.0,0.33333); norm /= 3.0; }
  Quark(const ParticleType& h,const QuantumProjections& qp,const Vektor3& p,const Vektor3& x)
    : RealParticle(h,qp,p,x,0) { reset(); L*=pow(3.0,0.33333); norm /= 3.0; }
  virtual ParticleBase* makeClone() const { return new Quark(*this); }
  virtual ~Quark();
  double getParam(int i) const { if (i==0) return path; else return justFiss; }
  void setParam(int i,double x) { i ? justFiss : path = x; }
  void refresh();
protected:
  virtual void announceEvent() const;
  virtual bool noPotentials() const { return false; }
  virtual double V(double r) const { return r; }
  virtual double dVdr(double r) const { return 1.0; }
};

class Hadron : public RealParticle
{
public:
  Hadron(const ParticleType& pp) : RealParticle(pp,-1) {}
  Hadron(const ParticleType& pp,const Vektor3& p,const Vektor3& x) 
    : RealParticle(pp,p,x,-1) {}
  Hadron(const ParticleType& pp,const QuantumProjections& qp ) 
    : RealParticle(pp,qp,-1) {}
  Hadron(const ParticleType& pp,const QuantumProjections& qp,const Vektor3& p,const Vektor3& x) 
    : RealParticle(pp,qp,p,x,-1) {}
  virtual bool isQuark() const { return false; }
  virtual bool isDiquark() const { return false; }
   
  virtual ParticleBase* makeClone() const { return new Hadron(*this); }
protected:
  virtual bool noPotentials() const { return true; }
  virtual double V(double r) const { return 0.0; }
  virtual double dVdr(double r) const { return 0.0; }
};

class ColorString 
{
  int N;
  double Etot;
  QuantumState pp;
  QuantumState array[2];
  Vektor3 beta,Ptot,Rtot;
public:
  ColorString(const double& Etot,const Vektor3& Ptot,const Vektor3& Rtot,int n,
	 const QuantumState parray[],const vector<ParticleBase*>&);
  ~ColorString();

  static double Cutoff;
  static double Ekin_min;
  void breakUp();
};

#include "Propagation.icc"

#endif

