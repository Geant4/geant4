#include <iostream.h>
#include <vector.h>
#include <string.h>
#include "math.hh"
#include "ParticleType.hh"
#include "reactionChannels.hh"
#include "String.hh"
#include "Random.hh"
#include "Memory.hh"
#include "DistributionFunction.hh"
#include "MathTools.hh"
#include "Collision.hh"

Knot<ParticleType>* Knot<ParticleType>::Root = new ParticleType;

RGB RGB::pickColor()
{ 
  return RGB(rand_gen::Random(1,3)); 
}

RGB operator+(const RGB& c1,const RGB& c2)
{
  int c;
  if ( c1*c2 > 0 && c1!=c2 ) {
    c = -sign((int)c1)*(abs(c1) ^ abs(c2));
  }
  else if ( c1 == 0 )
    c = c2;
  else if ( c2 == 0 )
    c = c1;
  else if ( c1 == -c2 )
    c = 0;
  else {
    throw "ERROR: wrong colors...";
  }
  return RGB(c);
}

RGB RGB::operator>>(int d)
{
  int x = abs(c)+d;
  while ( x>3 )
    x -= 3;
  RGB y = sign(c)*x;
  return y;
}

RGB RGB::operator<<(int d)
{
  int x = abs(c)-d;
  while ( x<1 )
    x += 3;
  RGB y = sign(c)*x;
  return y;
}

void QuantumNumbers::writeOut(ostream& o) const
{
  o << B() << "  " << S() << "  " << Isospin() << "  " << Spin() << "  ";
}

bool operator==(const ParticleType& x,const ParticleType& y) 
{
  return (x.id == y.id);
}

ParticleType::ParticleType() : minmass(0),width(0) {}

ParticleType::ParticleType(istream& in) 
  : minmass(0),width(0)
{
  static int id_counter = 0;
  String name;
  undef<double> u_mass,u_spin,u_iso,u_width,u_min;
  undef<int> u_nq,u_g,u_s,u_c;
  in >> name >> u_mass >> u_nq >> u_s >> u_c >> u_spin 
     >> u_iso >> u_g >> u_width >> u_min;
  g = u_g;
  if ( u_g.isValid() ) 
    setMaxColor((g>1) ? g : 1);
  if ( u_s.isValid() )
    setS(u_s);
  if ( u_c.isValid() )
    setC(u_c);
  if ( u_nq.isValid() && u_s.isValid() && u_c.isValid() )
    setB((u_nq - S()+Charm())/3.0);
  else if (u_nq.isValid())
    setB(u_nq/3.0);
  else
    setB(u_nq);
  if ( u_iso.isValid() )
    setIsospin(u_iso);
  if ( u_spin.isValid() )
    setSpin(u_spin);
  if ( u_mass.isValid() ) {
    setPeakMass(u_mass/1000.0); // GeV
    width = u_width/1000.0; // GeV
    BWnorm = width/(M_PI+2*atan(2*PeakMass()/width));
  }
  if ( u_min.isValid() ) 
    minmass = u_min;
  id = id_counter++;
  Insert(*this,name);
}

void ParticleType::getProbab(double m,vector<double>& P,vector<ParticleType*>& L, bool Integral) const
{
  if ( hasSuccessors() ) {
    for (int i=0; i<NSuccessors(); i++) 
      Successor(i).getProbab(m,P,L,Integral);
  }
  else {
    if ( m >= minmass || 1) {
    double x = P.size() ? P.back() : 0.0;
    double pp = 1;
    if ( width > 0 ) {
      if ( Integral ) 
	pp = 2.0/width*BWnorm*(atan(2*PeakMass()/width)+atan(2*(m-PeakMass())/width));
      else
	pp = BWnorm/(sqr(m-PeakMass())+sqr(width/2.0));
      P.insert(P.end(),x+pp);
      L.insert(L.end(),(ParticleType*)this);
    }
  }
  }
}

ParticleType::MassDist::MassDist(const ParticleType& h,double m) 
  : GaussIntegration(h.minmass),SimpleRejection(h.minmass,m),Type(h),m_max(m) 
{ 
  norm = 1.0/Integral(m); 
}

REAL ParticleType::MassDist::ToBeIntegrated(REAL x)
{
  vector<double> P;
  vector<ParticleType*> L;
  Type.getProbab(x,P,L);
  return P.back();
}

REAL ParticleType::MassDist::f(REAL x) const 
{ 
  return norm*((MassDist&)(*this)).ToBeIntegrated(x); 
}

double ParticleType::selectMass(double m_max) const
{
  MassDist Dist(*this,m_max);
  return Dist.getValue();
}

ParticleType& ParticleType::selectType(double m) const
{
  vector<double> P;
  vector<ParticleType*> L;
  getProbab(m,P,L,true);
  if ( P.empty() )
    throw "No suitable ParticleType found...";
  double r = rand_gen::Random()*P.back();
  int i=0;
  while ( P[i]<r ) 
    ++i;
  return *L[i];
}

ParticleType& ParticleType::selectType(int C_,const vector<ParticleBase*>& P,double m) const
{
  double y0 = 0;
  double i3= 0;
  double s3= 0;
  for (int j=0; j<P.size(); j++) {
    i3 += P[j]->Iso3();
    s3 += P[j]->Spin3();
  }
  vector<double> y;
  vector<int> list;
  vector<ParticleType*> L = getList();
  for (int i=0; i<L.size(); i++) {
    CollisionType& X = Knot<CollisionType>::FindKnot(*L[i]);
    double pc = CollisionType::FindDecomposition(C_,X,P);
    vector<double> P;
    vector<ParticleType*> Lnew;
    L[i]->getProbab(m,P,Lnew,false);
    double px = 0.0;
    if ( P.size() ) 
      px = P.back();
    if ( pc>0.0 && fabs(i3)<=L[i]->Isospin() && fabs(s3)<=L[i]->Spin() ) {
      y.insert(y.end(),y0+=pc*px);
      list.insert(list.end(),i);
    }
  }
  if ( y.size() ) {
    if ( y0 == 0.0 ) {
      cerr << "Too few energy: " << m << "!!\n";
      int i=0;
      while ( i<L.size() && L[i]->getWidth() == 0.0 ) { ++i; }
      if ( i==L.size() ) 
	return *L[0];
      else 
	return *L[i];
    }
    double r = rand_gen::Random()*y0;
    int l;
    for (l=0; l<y.size() && r>y[l]; l++);
    if (l==y.size()) {
      cerr << "Not found: 1\n";
//      printTree(cerr);
//      for (int i=0; i<P.size(); i++)
//	cerr << *P[i] << endl;
      throw "Not found";
    }
    return *L[list[l]];
  }
  else {
    cerr << "Undefined Particle!!\n";
    //    return Knot<ParticleType>::FindKnot("X");
//    printTree(cerr);
//    for (int i=0; i<P.size(); i++)
//      cerr << *P[i] << endl;
    throw undefinedParticle();  
  }
}

ParticleType& ParticleType::selectType() const
{
  vector<double> P;
  vector<ParticleType*> L = getList();
  int i = int(rand_gen::Random()*L.size());
  return *L[i];
}

double ParticleType::getMass(double Emax) const
{
  if ( joker() )
    throw "ERROR: getMass() called from Group...";
  double mass = PeakMass();
  if ( width>0.0 && Emax > 0.0 ) {
    mass = selectMass(Emax);
  }
  return mass;
}

double ParticleType::getLifetime() const
{
  if ( width > 0.0 ) {
    double x = -log(1.0-rand_gen::Random())/width;
    return x;
  }
  else
    return 1e+30;
}

double QuantumNumbers::getIso3() const 
{
  return int((2*Isospin()+1)*rand_gen::Random())-Isospin();
}

double QuantumNumbers::getSpin3() const 
{
  return int((2*Spin()+1)*rand_gen::Random())-Spin();
}

RGB QuantumNumbers::getColor(RGB col) const 
{
  int c = col;
  if ( !c ) 
    c = (maxColor()>1) ? int(maxColor()*rand_gen::Random())+1 : 0;
  if ( c && fabs(B()) > 0.5 )
    c = -c;
  return RGB(c);
}

void ParticleType::writeOut(ostream& o) const
{
  o << Name() << ": " << PeakMass() << "  ";
  QuantumNumbers::writeOut(o);
}

double ParticleType::isEqualTo(const ParticleType& x) const 
{ 
  double succ = matches(x);
  return succ;
}

