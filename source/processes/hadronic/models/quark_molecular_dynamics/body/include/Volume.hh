#ifndef __VOLUME__
#define __VOLUME__

#include "ThermDist.hh"
#include "ParticleType.hh"
#include "HadronGas.hh"

class Particle;
class Geometry;

template<class t>
class Volume 
{
  typedef t Type;
  REAL T,mu,mu_s;
  grandCanonicalEnsemble Ensemble;
protected:
  vector<ParticleType*>& List;
  Geometry& G;
public:
  Volume(Geometry& G_,vector<ParticleType*>& L,REAL T_,REAL mu_,REAL mus_,int num = 0);
  double Etot() const { return Ensemble.Etot(parameter(T,G.getVolume(),mu,mu_s)); }
  void createParticles(double frac = 0.5);
private:
  virtual void setSpecies(int i,int n,int Nbar,double frac = 0.5);
};

template<class t>
class QuarkVolume : public Volume<t>
{
public:
  QuarkVolume(Geometry& G_,vector<ParticleType*>& L,REAL T_,REAL mu_,REAL mus_,int num = 0)
    : Volume<t>(G_,L,T_,mu_,mus_,num) {}
private:
  void setSpecies(int i,int n,int Nbar,double frac = 0.5);
};


#include "Volume.tcc"

#endif
