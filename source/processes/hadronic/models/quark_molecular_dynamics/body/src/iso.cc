#include "g4std/vector"
#include "math.hh"
#include "iso.hh"
#include "clebsch.hh"
#include "Random.hh"
#include "MultiIndex.hh"
#include "array.hh"
#include "ParticleKinematics.hh"

void Iso::Projections(int n,double J,double M,double* jk,double* mk,bool* isSet)
{
  vector<MultiIndex*> all;
  vector<double> probab;
  Array<int> jkmax(n);
  {for ( int i=0; i<n; i++)
    mk[i] = 0;}
  if ( !isSet ) {
    isSet = new bool[n];
    for (int i=0; i<n; i++)
      isSet[i] = false;
  }
  for (int i=0; i<n; i++)
    if ( !isSet[i] )
      jkmax[i] = int(2*jk[i]+1);
    else 
      jkmax[i] = 1;
  int num = 0;
  double Ptot = 0.0;
  for (MultiIndex proj(n,jkmax); proj.isValid(); proj++) {
    double mtot = 0.0;
    for ( int k=0; k<n; k++)
      if ( !isSet[k] ) 
	mtot += proj[k]-jk[k];
      else
	mtot += mk[k];
    if ( mtot == M ) {
      all.insert(all.end(),new MultiIndex(proj));
      for (int i=0; i<n; i++)
	if ( !isSet[i] ) 
	  mk[i] = -jk[i]+proj[i];
      double P = sqr(f(n-2,jk[n-1],mk[n-1],J,M,jk,mk));
      //      G4cerr << proj << "  " << P << G4endl;
      Ptot += P;
      probab.insert(probab.end(),P);
      ++num;
    }
  }
  if ( probab.empty() )
    return;
  double r = Ptot*rand_gen::Random();
  double s = probab[0];
  int j = 0;
  for (j=0; j<num && r>s; j++) 
    s += probab[j+1];
  for (int k=0; k<n; k++)
    if ( !isSet[k] )
      mk[k] = (*all[j])[k]-jk[k];
}

double Iso::f(int n,double jnp1,double mnp1,double J,double M,double* jk,double* mk) {
  double s = 0;
  if ( n ) {
    double jmin = fabs(jk[n]-jnp1);
    double jmax = jk[n]+jnp1;
    double mn = mk[n]+mnp1;
    for ( double jn=jmin; jn<=jmax; jn++) {
      s += ClebschGordan::Clebsch(jk[n],mk[n],jnp1,mnp1,jn,mn)*f(n-1,jn,mn,J,M,jk,mk);
    }
  }
  else
    s = ClebschGordan::Clebsch(jk[n],mk[n],jnp1,mnp1,J,M);
  return s;
}

double Iso::chooseMultiplett(double j1,double m1,double j2,double m2,states which)
{
  int n = int(j1+j2-fabs(j1-j2)+1);
  Array<double> probJ(n);
  double x = 0.0;
  double J = fabs(j1-j2);
  for (int i=0; i<n; i++) {
    x += sqr(ClebschGordan::Clebsch(j1,m1,j2,m2,J,m1+m2));
    if ( which == LOWEST && x>0 )
      return J;
    probJ[i] = x;
    J += 1.0;
  }
  double r = rand_gen::Random()*x;;
  int k;
  for (k=0; k<n && r>probJ[k]; k++);
  return fabs(j1-j2)+k;
}

