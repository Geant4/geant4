#include "ParticleKinematics.hh"
#include "ParticleProperties.hh"
#include "G4ios.hh"
#include "math.hh"
#include "Random.hh"
#include "iso.hh"
#include "Collision.hh"
#include "EventHandling.hh"


QuantumProjections::QuantumProjections(int c_,RGB col,double i3,double s3)
  : c(c_),color(col),iso3(i3),spin3(s3) {}

QuantumProjections::QuantumProjections(const ParticleType& h)
  : c(1),color(h.getColor()),iso3(h.getIso3()),spin3(h.getSpin3()) {}

QuantumProjections& QuantumProjections::operator=(const QuantumProjections& Q)
{
  c = Q.c; color = Q.color; iso3 = Q.iso3; spin3 = Q.spin3; 
  return *this;
}

QuantumState& QuantumState::operator=(const QuantumState& x)
{
  QuantumProjections::operator=(x);
  Pattern<double>::operator=(x);
  return *this;
}

QuantumState anti(const QuantumState& p0) 
{
  QuantumState p = p0;
  p.c = -p0.c;
  p.color = -p0.color;
  p.iso3 = -p0.iso3;
  p.spin3 = -p0.spin3;
  return p;
}

QuantumState operator+(const QuantumState& p1,const QuantumState& p2)
{
  QuantumState p;
  if ( p1.C() == p2.C() ) 
    p.c = p1.C();
  else
    p.c = 1;
  p.color = p1.color+p2.color;
  p.iso3 = p1.iso3+p2.iso3;
  p.spin3 = p1.spin3+p2.spin3;

  if ( fabs(p1.iso3+p2.iso3) == p1.Isospin()+p2.Isospin() ) 
    p.setIsospin(p1.Isospin()+p2.Isospin());
  if ( fabs(p1.spin3+p2.spin3) == p1.Spin()+p2.Spin() ) 
    p.setSpin(p1.Spin()+p2.Spin());
  p.setB(p.C()*(p1.C()*p1.B()+p2.C()*p2.B()));
  p.setS(p.C()*(p1.C()*p1.S()+p2.C()*p2.S()));
  p.setC(p.C()*(p1.C()*p1.Charm()+p2.C()*p2.Charm()));
  p.setMaxColor(p.color ? 3 : 1);
  return p;
}


QuantumState addToLowest(const QuantumState& p1,const QuantumState& p2)
{
  QuantumState p;
  if ( p1.C() == p2.C() ) 
    p.c = p1.C();
  else
    p.c = 1;
  p.color = p1.color+p2.color;
  p.iso3 = p1.iso3+p2.iso3;
  p.spin3 = p1.spin3+p2.spin3;

  p.setIsospin(Iso::chooseMultiplett(p1.Isospin(),p1.iso3,p2.Isospin(),p2.iso3,Iso::LOWEST));
  p.setSpin(Iso::chooseMultiplett(p1.Spin(),p1.spin3,p2.Spin(),p2.spin3,Iso::LOWEST));
  p.setB(p.C()*(p1.C()*p1.B()+p2.C()*p2.B()));
  p.setS(p.C()*(p1.C()*p1.S()+p2.C()*p2.S()));
  p.setC(p.C()*(p1.C()*p1.Charm()+p2.C()*p2.Charm()));
  p.setMaxColor(p.color ? 3 : 1);
  return p;
}

ParticleProperties::ParticleProperties(const ParticleType& h,double Emax)
  : QuantumProjections(1,h.getColor(),h.getIso3(),h.getSpin3()),Type(h),
    mass(h.getMass(Emax)),flag(0)
{}

ParticleProperties::ParticleProperties(const ParticleType& h,const QuantumProjections& q,double Emax)
  : QuantumProjections(q),Type(h),
    mass(h.getMass(Emax)),flag(0)
{
}

ParticleProperties::ParticleProperties(const ParticleProperties& p)
  : QuantumProjections(p),Type(p.Type),
    mass(p.mass),lifetime(p.lifetime),flag(p.flag)
{
}

void ParticleProperties::writeOut(G4std::ostream& o) const
{
  Type.writeOut(o);
  o << c << "  " << color << "  " << iso3 << "  " << spin3 << "  ";
}

void ParticleProperties::ChargeProjection(double ch)
{ 
  iso3 = ch-(Type.B()+Type.S())/2.0; 
}

String ParticleProperties::Name() const 
{ 
  String s;
  if ( C()<0 && B() ) 
    s = "anti-";
  s += Type.Name()+printCharge();
  return s;
}

String ParticleProperties::printCharge() const
{
  String s;
  //  int ch = int((B() ? sign(B()) : 1)*Charge());
  int ch = (int)Charge();
  char c;
  if (ch>0) c = '+';
  if (ch<0) c = '-';
  if (ch==0) c = '0';
  int j=0;
  do {
    s+=c;
    ++j;
  }
  while (j<abs(ch));
  if ( j>2 ) {
    int ch1 = (int)Charge();
    G4cerr << "ATTENTION: " << ch << "  " << c << G4endl;
  }
  return s;
}

void ParticleProperties::setLifetime(double lt) 
{ 
  if ( lt < 0 ) 
    lifetime = Time()+Type.getLifetime(); 
  else if ( Type.getWidth() > 0 ) {
    vector<CollisionTab*>::iterator Y = CollisionTab::exists(this);
    if ( Y != CollisionTab::Root.end() ) 
      CollisionTab::remove(Y);
    lifetime = lt;
  }
  else 
    lifetime = 1e30;
  if ( lifetime<1e30 ) {
    CollisionTab::addEntry(lifetime,*this,ALL); 
  }
}
