#include "G4ios.hh"
#include "metropolis.hh"
#include "Random.hh"
#include "Geometry.hh"

G4std::ostream& operator<<(G4std::ostream& o,const Metropolis& M) {
  for (int i=0; i<M.N; i++) {
    o << M.getX(i)[1] << "  " << M.getX(i)[2] << "  " << M.getX(i)[3] << "  "
      << M.getP(i)[1] << "  " << M.getP(i)[2] << "  " << M.getP(i)[3] << G4endl;
  }
  return o;
}

Metropolis::Metropolis(Geometry& G_,int n,double T_,double a_,double b_) 
  : N(n),N_init(-1),T(T_),E_k(0.0),F(0.0), 
    k(0),kstart(1000),l(0),a(a_),b(b_),G(G_)
{
}

Metropolis::~Metropolis()
{
  //  delete [] x;
  //  delete [] p;
}

double Metropolis::E_n(int n) const
{
  double e = E_part(n);
  for (int i=0; i<N; i++)
    e += E_int(n,i);
  return e;
}

void Metropolis::evaluate(long MAX_ITER) 
{
  int M = 300,m0=0;
  double F = 0, F2 = 0;
  double S = 0,oldS=0,dS=0;
  G4cerr << N << G4endl;
    long count = 0;
    k = 0;
    //    while ( fabs(S-oldS)>=dS ) {
    while ( ++count<MAX_ITER ) {
      ++k;
      if ( k>kstart ) {
	/*
	if ( k % M ) {
	  ++m0;
	  double f = Observable();
	  F += f;
	  F2 += f*f;
	}
	else {
	  oldS = S;
	  S = F/m0; 
	  dS = sqrt((F2-m0*S*S)/(m0-1));
	  //	  G4cerr << k << "  " << S << "  " << oldS << "  " << fabs(S-oldS) << "  " << dS << G4endl;
	  F = 0;
	  F2 = 0;
	}
	*/
      }
      if ( !(k % N) ) check();
      if ( k>kstart && !(k % N) ) {
	doSomething();
      }
      int n = int(N*rand_gen::Random());
      Vektor3 x0 = getX(n);
      Vektor3 p0 = getP(n);
      double E_old = E_n(n);
      Vektor3 xx;
      double xi;
      do {
	 xi = rand_gen::Random();
	xx = getX(n)+Vektor3::isotropy(xi*a);
      }
      while ( ! G.isInside(xx) );
      setX(n,xx);
      setP(n,p0+Vektor3::isotropy(xi*b));
      double dE = E_n(n)-E_old;
      double w = exp(-dE/T);
      if ( rand_gen::Random() > w ) {
	setX(n,x0);
	setP(n,p0);
      }
      else {
	E_old += dE;
      }
      /*
      {
	Vektor3 p0 = getP(n);
	double xi = rand_gen::Random();
	setP(n,p0+Vektor3::isotropy(xi*b));
	double dE = E_n(n)-E_old;
	double w = exp(-dE/T);
	if ( rand_gen::Random() > w ) {
	  setP(n,p0);
	}
	else {
	  //      E_old += dE;
	}
      }
      */
    }
}

void Metropolis::setParticle(const Vektor3& x_,const Vektor3& p_)
{
  if ( N_init<N-1 ) {
    ++N_init;
    setX(N_init,x_);
    setP(N_init,p_);
  }
}

