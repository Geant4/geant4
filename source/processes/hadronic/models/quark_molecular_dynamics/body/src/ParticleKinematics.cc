#include "math.hh"
#include "newvector.hh"
#include "ParticleKinematics.hh"

double CMmomentum(double M,double m1,double m2)
{
  if ( M>=m1+m2 ) 
    return sqrt(sqr(M*M-m1*m1-m2*m2)
		-4.0*sqr(m1*m2))/(2.0*M);
  else
    throw "M<m1+m2";
}

Vektor4 Lorentz(const Vektor3& beta,const Vektor4& P)
{
  double b2 = beta*beta;
  if ( b2>0.0 ) {
    double gamma=1.0/sqrt(1.0-b2);
    Vektor3 Px = spatial(P);
    return Vektor4(Px + ((gamma-1)/b2*(beta*Px))*beta-(gamma*P[0])*beta,gamma*(P[0]-beta*Px));
  }
  else
    return P;
}
