#include "HadronGas.hh"
#include "Random.hh"


const REAL gasEnsemble::volume = 1.0/pow(0.197,3);
REAL RelBoltzmann::norm = 0.5/sqr(3.1415926)*gasEnsemble::volume;

void parameter::init(int n_) 
{ 
  n = n_;
  N.reset(2*n_);
  for (int i=1; i<=dim(N); i++) {
    N[i] = 0.0;
  }
  pressure = 0;
}

void parameter::init(int n_,const Vektor& N_) 
{ 
  n = n_;
  N = N_;
}

parameter& parameter::operator=(const parameter& p)
{
  T = p.T;
  mu_q = p.mu_q;
  mu_s = p.mu_s;
  V = p.V;
  if ( p.n ) {
    N = p.N;
  }
  else
    N.reset(0);
  return *this;
}

parameter& parameter::operator+=(const parameter& p)
{
  if ( p.n ) {
    N+= p.N;
  }
  else {
    throw "parameter::operator+=";
  }
  return *this;
}

parameter parameter::operator*(REAL y)
{
  parameter p(*this);
  if ( p.n ) {
    N *= y;
  }
  else {
    throw "parameter::operator*";
  }
  return p;
}


REAL gasEnsemble::SumN(const parameter& p,signed int B) const
{
  REAL S=0.0;
  for (vector<ParticleType*>::iterator X=List.begin(); X!=List.end(); X++) {
    S += N(*(*X),p,B);
  }
  return S;
}

Vektor grandCanonicalEnsemble::q_tot(const parameter& p) const
{
  Vektor vars(3);
  for (vector<ParticleType*>::iterator X=List.begin(); X!=List.end(); X++) {
    REAL lambda = K->TdlnZ((*X)->getMass(),p.T);
    REAL netn = NetN(*(*X),p);
    REAL n = N(*(*X),p);
    vars[1] += netn*(*X)->Nq();
    vars[2] += netn*(*X)->S();
    vars[3] += p.T*n*lambda;
  }
  return vars;
}

REAL grandCanonicalEnsemble::Mu(const ParticleType& h,const parameter& p,signed int B) const
{ 
  REAL y,shift = p.T*log(h.degeneracy());
  if ( B==0 )
    y = 2*(h.Nq()*p.mu_q+h.S()*p.mu_s);
  else {
    y = B*(h.Nq()*p.mu_q+h.S()*p.mu_s);
  }
  return y;
}

REAL gasEnsemble::Etot(const parameter& p) const
{
  REAL e = 0.0,zeta;
  for (vector<ParticleType*>::iterator X=List.begin(); X!=List.end(); X++) {
    zeta = K->TdlnZ((*X)->getMass(),p.T);
    e += N(*(*X),p,0)*zeta;
  }
  return p.T*e;
}

REAL grandCanonicalEnsemble::N(const ParticleType& h,const parameter& p,signed int B) const
{
  int xx = h.Nq();
  REAL mu = p.mu_s*h.S() + p.mu_q*h.Nq();
  REAL z;
  //  if ( h.B != 0.0 )
  //    B = abs(B);
  if ( B==0 ) 
    z = exp((mu)/p.T)+exp((-mu)/p.T);
  else
    z = exp((B*mu)/p.T);
  return p.V*h.degeneracy()*z*K->Z(h.getMass(),p.T);
}

REAL  gasEnsemble::Entropy(const parameter& p) const
{
  REAL S = 0.0,zeta;
  for (vector<ParticleType*>::iterator X=List.begin(); X!=List.end(); X++) {
    zeta = K->TdlnZ((*X)->getMass(),p.T);
    S += N(*(*X),p,0)*(1+zeta);
  }
  return S;
}

REAL canonicalEnsemble::N(const ParticleType& h,const parameter& p,signed int B)
{
  int i;
	int ListSize = List.size();
  for (i=0; i<ListSize && *List[i] != h; i++);
  switch ( B ) {
  case -1 : return p.N[i+ListSize]; break;
  case 1 : return p.N[i]; break;
  case 0 : return p.N[i]+p.N[i+ListSize]; 
  }
}

parameter&  grandCanonicalEnsemble::initGas(parameter& p)
{
	int i;
	int ListSize = List.size();
  p.init(ListSize);
  for (i=0; i<ListSize; i++) {
    p.N[i+1] = N(*List[i],p,1);
    p.N[i+ListSize+1] = N(*List[i],p,-1);
  }
  p.pressure = SumN(p)*p.T/p.V;
  return p;
}
