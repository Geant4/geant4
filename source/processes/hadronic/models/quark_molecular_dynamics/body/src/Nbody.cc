#include "Nbody.hh"
#include "reactionChannels.hh"
#include "ParticleType.hh"
#include "Collision.hh"

Nbody* Particle::Soup = 0;
double Particle::rho_0 = 0.16; // fm^{-3}

Nbody::Nbody(double h_) 
  : time(0.0),h(h_),array(0),forces(0),firstCall(true),N(0),Npart(0) {}


Nbody::~Nbody() 
{
  vector<Particle*>::iterator X=List.begin(); 
  while ( !List.empty() ) {
    delete List.back();
  }
}

void Nbody::add(Particle* p,int i)
{
  int where = N;
  int k = List.size();
  if ( i>=0 && k>0 ) {
    where = List[i]->offset;
    k = i;
  }
  for (unsigned int j=k; j<List.size(); j++)
    List[j]->offset += 6;
  List.insert(List.begin()+k,p);
  array.reset(6,where);
  forces.reset(6,where);
  p->offset = where;
  p->time = time;
  N += 6;
  ++Npart;
}

void Nbody::sub(Particle* p)
{
  int i=0;
  vector<Particle*>::iterator X = List.begin();
  while ( *X != p && X != List.end() ) { ++i; ++X; }
  if ( X!=List.end() ) {
    for (vector<Particle*>::iterator Y = X+1; Y!=List.end(); Y++) 
      (*Y)->offset-=6;
    sub(i);
  }
  List.erase(X);
}

void Nbody::sub(int i)
{
  int pos = List[i]->offset;
  array.remove(pos+1,6);
  forces.remove(pos+1,6);
  N -= 6;
  --Npart;
}


void Nbody::handleCollisions()
{
  /*
  for (int i=0; i<Npart; i++)
    for (int j=i+1; j<Npart; j++) {
      vector<ParticleBase*> L;
      L.insert(L.end(),(ParticleBase*)List[i]);
      L.insert(L.end(),(ParticleBase*)List[j]);
      CollisionType* C = CollisionType::checkCollision(L);
      if ( C ) 
	   CollisionTab::addEntry(time,L);
    }
    */
  //    List[i]->announceEvent();
  CollisionTab::perform(time);
}

Vektor& Nbody::function(Vektor& f)
{
  for ( int i=0; i<Npart; i++ ) {
    Vektor3 dx = dHdp(i);
    Vektor3 dp = dHdx(i);
    for (int k=1; k<=3; k++) {
      f[List[i]->offset+k] = dx[k];
      f[List[i]->offset+k+3] = -dp[k];
    }
    List[i]->rho = 0;
    for (int j=0; j<Npart; j++) 
      List[i]->rho += List[i]->gauss(List[j]->Coordinates());
  }
  return f;
}


void Nbody::one_step() 
{
  if ( h>0 ) 
  {
    Vektor f0(array);
    Vektor k1(N),k2(N),k3(N);
    
    if ( firstCall ) {
      function(forces);
      firstCall = false;
    }
    array += (0.5*h)*forces;
    function(k1);
    array = f0+(0.5*h)*k1;
    function(k2);
    array = f0+h*k2;
    function(k3);
    time += h;
    for (int k=1; k<=N; k++) {
      array[k] = f0[k]+(h/6.0)*(forces[k]+2.0*k1[k]+2.0*k2[k]+k3[k]);
    }
  }
  checkRange();
  function(forces);
  handleCollisions();
}


void Nbody::print(G4std::ostream& o) 
{ 
  o << "#  time=" << Time() << ", Npart=" << Npart << G4endl;
  for (int i=0; i<Npart; i++) { o << *List[i] << G4endl; }
}


Particle::Particle() : L(0.287),norm(rho_0)
{ 
  Soup->add(this,-1); 
}

Particle::Particle(int pos) : L(0.287),norm(rho_0)
{
  Soup->add(this,pos);
}

Particle::Particle(const Vektor3& mom,const Vektor3& x,int pos) 
  : L(0.287),norm(rho_0)
{
  Soup->add(this,pos);
  SetCoordinates(x);
  SetMomentum(mom);
}

Particle::~Particle()
{
  Soup->sub(this);
}


Vektor3 Particle::Vprime(const ParticleBase& y)
{
  Vektor3 d = y.Coordinates()-Coordinates();
  double ld = length(d);
  return (dVdr(ld)/ld)*d;
}

double Particle::gauss(const Vektor3& x) const
{
  return norm*exp(-0.5*sqr(length(Soup->dr(x,Coordinates()))/L));
}
