#ifdef IS_GCC
#pragma implementation "Geometry.hh"
#endif

#include "Random.hh"
#include "Geometry.hh"

bool Geometry::moveToSurface(Vektor4& x,const Vektor4& p) {
  bool ret = true;
  Vektor3 p0 = spatial(p);
  Vektor3 velocity = p0/p[0];
  Vektor3 x0 = spatial(x);
  Vektor3 r0;
  try {
    REAL t = hitSurface(velocity,x0,r0);
    x[0] += t;
    x.setSpatial(r0);
  }
  catch ( ... ) { ret = false; }
  return ret;
}

bool Geometry::moveToSurface(Vektor4& x,Vektor4& p,Matrize& U,Vektor3& r0) {
  Vektor3 p0 = spatial(p);
  Vektor3 velocity = p0/p[0];
  Vektor3 x0 = spatial(x);
  try {
    REAL t = hitSurface(velocity,x0,r0);
    Vektor3 null;
    U = rotateFrame(r0);
    p.setSpatial(U*p0);
    x[0] += t;
    x.setSpatial(null);
    return true;
  }
  catch (...) { return false; }
}


REAL halfSpace::frontArea = 36; // fm^2

REAL halfSpace::hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0)
{
  if ( v[3] <= 0 ) 
    throw "Never reaches surface...";
  REAL tmin = -x0[3]/v[3];
  r0 = x0+v*tmin;
  for (int j=1; j<=2; j++)
    while ( fabs(r0[j]) > 0.5*sqrt(frontArea) )
      r0[j] -= sign(r0[j])*sqrt(frontArea);
  return tmin;
}

bool halfSpace::homogeneous(Vektor4& x)
{
  REAL l0 = (dL != 0.0) ? dL : L;
  x[0] = 0.0;
  for (int i=1; i<=2; i++)
    x[i] = (rand_gen::Random()-0.5)*sqrt(frontArea);
  x[3] = -rand_gen::Random()*l0;
  return true;
}

Matrize halfSpace::rotateFrame(const Vektor3&)
{
  Matrize U(3,3);
  for (int i=1; i<=3; i++)
    U(i,i) = 1.0;
  return U;
}

Matrize Box::rotateFrame(const Vektor3& r)
{
  Matrize U(3,3);
  signed int j=0;
  for (int i=1; i<=3 && j==0; i++) {
    //    G4cerr << i << "  " << fabs(r[i]/(0.5*L))-1.0 << "  " << (signed int)(r[i]/(0.5*L)) << G4endl;
    j=round(r[i]/(0.5*L))*i;
  }
  switch ( j ) {
  case 3 : U(1,1) = 1; U(2,2) = 1; U(3,3) = 1; break;
  case -3 : U(1,1) = 1; U(2,2) = -1; U(3,3) = -1; break;
  case 1 : U(1,3) = -1; U(2,2) = 1; U(3,1) = 1; break;
  case -1 : U(1,3) = 1; U(2,2) = 1; U(3,1) = -1; break;
  case 2 : U(1,1) = 1; U(2,3) = 1; U(3,2) = 1; break;
  case -2 : U(1,1) = 1; U(2,3) = -1; U(3,2) = -1; break;
  default: G4cerr << "Fehler!!! " << j << G4endl; exit(0);
  }
  return U;
}


void Box::reflect(Vektor3& x) const
{
  REAL d = 0.5*L;
  for (int k=1; k<=3; k++) {
    if ( x[k]>d ) {
      x[k] -= L;
    }
    else if ( x[k] < -d ) {
      x[k] += L;
    }
  }
}

REAL Box::hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0)
{
  REAL tmin = 1e10;
  for (int i=1; i<=3; i++) {
    REAL t = (0.5*L*sign(v[i])-x0[i])/v[i];
    if ( t<tmin) 
      tmin = t;
  }
  r0 = x0+v*tmin;
  return tmin;
}

bool Box::homogeneous(Vektor4& x)
{
  REAL l0 = (dL != 0.0) ? L-dL : REAL(0.0);
  x[0] = 0.0;
  for (int i=1; i<=3; i++) {
    REAL y = (rand_gen::Random()-0.5);
    x[i] = sign(y)*(l0+fabs(y)*(L-l0));
  }
  return true;
}

Vektor3 InfiniteBox::dr(const Vektor3& x1,const Vektor3& x2) const
{
  Vektor3 d = x1-x2;
  for (int k=1; k<=3; k++)
    if ( fabs(d[k]) > 0.5*size() ) 
      d[k] -= sign(d[k])*size();
  return d;
}

bool Sphere::homogeneous(Vektor4& x)
{
  REAL l0 = (dR != 0.0) ? R-dR : REAL(0.0);
  REAL r = l0+(R-l0)*pow(rand_gen::Random(),1.0/3.0);
  x.setSpatial(Vektor3::isotropy(r));
  return true;
}


Matrize Sphere::rotateFrame(const Vektor3& r)
{
  Matrize U(3,3);
  Vektor3 r0 = 1.0/sqrt(square(r))*r;
  REAL rho = sqrt(sqr(r0[1])+sqr(r0[2]));
  signed int s = 1;
  
  for (int i=1; i<=3; i++) 
    U(3,i) = r0[i];
  U(2,1) = r0[2]/rho*s;
  U(2,2) = -r0[1]/rho*s;
  U(1,1) = -r0[1]*r0[3]/rho*s;
  U(1,2) = -r0[2]*r0[3]/rho*s;
  U(1,3) = rho*s;
  /*
  REAL r2 = square(r);
  REAL r1 = sqrt(r2);
  REAL cos_theta = r[3]/r1;
  REAL sin_theta = sqrt(1-sqr(cos_theta));
  REAL sin_phi = r[2]/(r1*sin_theta);
  REAL cos_phi = r[1]/(r1*sin_theta);
  
  U(1,1) = cos_phi*cos_theta;
  U(1,2) = sin_phi*cos_theta;
  U(1,3) = sin_theta;
  U(2,1) = -sin_phi;
  U(2,2) = cos_phi;
  U(3,1) = -cos_phi*sin_theta;
  U(3,2) = -sin_phi*sin_theta;
  U(3,3) = cos_theta;
  */
  
  return U;
}

REAL Sphere::hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0)
{
  REAL v2 = square(v);
  REAL x2 = square(x0);
  REAL cos_chi = (v*x0)/sqrt(v2*x2);
  REAL sin_chi = 1-sqr(cos_chi);
  REAL tmin = (sqrt(R*R-x2*sin_chi) - sqrt(x2)*cos_chi)/sqrt(v2);

  r0 = x0+v*tmin;
  return tmin;
}

bool Tube::homogeneous(Vektor4& x)
{
  REAL r = R*pow(rand_gen::Random(),0.5);
  REAL Phi = 2*mathConstants::Pi*rand_gen::Random();
  x[1] = r*cos(Phi);
  x[2] = r*sin(Phi);
  x[3] = (rand_gen::Random()-0.5)*L;

  return true;
}


Matrize Tube::rotateFrame(const Vektor3& r)
{
  Matrize U(3,3);
  Vektor3 r0;
  REAL rho2 = sqr(r[1])+sqr(r[2]);
  REAL rho = sqrt(rho2);
  if ( fabs(r[3]) < 0.5*L ) {
    r0[1] = R/rho*r[1];
    r0[2] = R/rho*r[2];
    r0[3] = r[3];
    U(1,3) = -1;
    U(2,1) = -r0[2]/R;
    U(2,2) = r0[1]/R;
    U(3,1) = r0[1]/R;
    U(3,2) = r0[2]/R;
  }
  else {
    REAL dz = r[3]-sign(r[3])*0.5*L;
    REAL u = sqrt(rho2+sqr(dz));
    r0[1] = 1/u*r[1];
    r0[2] = 1/u*r[2];
    r0[3] = 1/u*dz;
    rho = sqrt(sqr(r0[1])+sqr(r0[2]));
    signed int s = 1;
    
    for (int i=1; i<=3; i++) 
      U(3,i) = r0[i];
    U(2,1) = r0[2]/rho*s;
    U(2,2) = -r0[1]/rho*s;
    U(1,1) = -r0[1]*r0[3]/rho*s;
    U(1,2) = -r0[2]*r0[3]/rho*s;
    U(1,3) = rho*s;
  }
  
  return U;
}

REAL Tube::hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0)
{
  REAL v2 = sqr(v[1])+sqr(v[2]);
  REAL x2 = sqr(x0[1])+sqr(x0[2]);
  REAL vx = x0[1]*v[1]+x0[2]*v[2];
  REAL cos_chi = (vx)/sqrt(v2*x2);
  REAL sin_chi = 1-sqr(cos_chi);
  REAL tmin = (sqrt(R*R-x2*sin_chi) - sqrt(x2)*cos_chi)/sqrt(v2);
  if ( fabs(x0[3]+v[3]*tmin) > L/2 ) {
    REAL t1 = 0;
    if ( fabs(x0[3]) < 0.5*L)
      t1 = (0.5*L*sign(v[3])-x0[3])/v[3];
    Vektor3 x1 = x0+v*t1;
    x1[3]-=sign(x1[3])*0.5*L;
    REAL x12 = square(x1);
    v2 = square(v);
    cos_chi = (v*x1)/sqrt(v2*x12);
    sin_chi = 1-sqr(cos_chi);
    REAL t2 = (sqrt(R*R-x12*sin_chi) - sqrt(x12)*cos_chi)/sqrt(v2);
    tmin = t1+t2;
  }
  r0 = x0+v*tmin;
  return tmin;
}




