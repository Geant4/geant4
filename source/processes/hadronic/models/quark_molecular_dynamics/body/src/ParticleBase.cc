#include "ParticleBase.hh"
#include "String.hh"
#include "ParticleType.hh"
#include "math.hh"

G4std::ostream& operator<<(G4std::ostream& o,const ParticleBase& p)
{
  o << p.Name()        << " (" 
    << p.PDGCode()     << "): " 
    << p.B()           << "  " 
    << p.S()           << "  "
    << p.Charm()       << "  " 
    << p.Mass()        << "  " 
    << p.Isospin()     << "  " 
    << p.Iso3()        << "  " 
    << p.Spin()        << "  " 
    << p.Spin3()       << "  " 
    << p.Coordinates() << "  " 
    << p.Momentum()    << "  " 
    << p.Color()       << "  " 
    << p.Force()       << "  " 
    << p.Lifetime()    << "  " 
    << p.Flag()        << "  " ;
  return o;
}

void ParticleBase::Lorentz(int j,double beta)
{
  double gamma=1.0/sqrt(1.0-beta*beta);
  double pj = Momentum(j);
  Momentum(j) = gamma*(pj-beta*E());
}

void ParticleBase::Lorentz(const Vektor3& beta)
{
  double b2 = beta*beta;
  double gamma=1.0/sqrt(1.0-b2);
  SetMomentum(Momentum()+ ((gamma-1)/b2*(beta*Momentum()))*beta-(gamma*E())*beta);
}

double ParticleBase::TransverseMass() const
{ 
  return sqrt(sqr(Mass())+sqr(Momentum(1))+sqr(Momentum(2))); 
}

double ParticleBase::TransverseMomentum() const
{ 
  return sqrt(sqr(Momentum(1))+sqr(Momentum(2))); 
}
