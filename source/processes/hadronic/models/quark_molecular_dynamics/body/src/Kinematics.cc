#include <iomanip.h>
#include <algo.h>
#include "Kinematics.hh"

REAL FragmentationFunction::a = 1.0;
REAL FragmentationFunction::b = 0.7; // Units: [GeV^(-2)]

REAL Banerjee::zeta = 0.0;

Kinematics::chooseKin Kinematics::whichKin = &Kinematics::kin1;

void Kinematics::calculate(double sigma)
{
  p1 = p0; p2 = p0;
  x1 = x0; x2 = x0;
  sigma_l = sigma*p0[3];

  REAL phi = 2*mathConstants::Pi*rand_gen::Random();
  // Determine transverse motion
  GaussDeviates PT;
  REAL Qt = 0.0;
  int counter = 100;
  do {
    Qt = sqrt(sigma/2.0/mathConstants::Pi)*fabs(PT.getValue());
  }
  while ( Qt > sigma_l && --counter ); 
  if ( !counter ) {
    throw "Counter overflow...";
    //    throw NotPossible(p1->Id());
  }
  p1[1] = Qt*sin(phi);
  p1[2] = Qt*cos(phi);
  p2[1] = p0[1]-p1[1];
  p2[2] = p0[2]-p1[2];

  pq2 = sqr(p0[3]);
  pq = sqrt(pq2);
  REAL Qt2 = sqr(Qt);
  REAL Pt2 = sqr(p2[1])+sqr(p2[2]); // new
  Mt2 = (sqr(m2) + Pt2)/pq2;
  Mt = sqrt(Mt2);
  mQt2 = (sqr(m1)+Qt2)/pq2;
  mQt = sqrt(mQt2);

  eps = (p0[0]+sigma_l)/pq;
  eps_min = sqrt(1+sqr(Mt+mQt));
  eps_max_0 = Mt+sqrt(1+mQt2);
  eps_max_1 = mQt+sqrt(1+Mt2);

  if ( eps < eps_min ) {
    throw "Not possible...";
  }

  REAL z=0.0, l=0.0;

  ((*this).*whichKin)(z,l);


  if ( l<0 || l>sigma_l ) {
    cerr << "Fehler: " << z << "  " << l << "  " 
	 << pq*eps_min-p0[0] << "  " 
	 << sigma_l << endl;
    throw "Not possible...";
  }

  p1[3] = (1-z)*p0[3];
  p1[0] = sqrt(m1*m1+sqrt(p1[3]));
  p2[3] = z*p0[3];
  p2[0] = sqrt(m2*m2+sqrt(p2[3]));
  x1[3] -= l/sigma;
}

REAL Kinematics::get_eps(REAL z) {
  return pq*(sqrt(z*z+Mt2)+sqrt(sqr(1-z)+mQt2))-p0[0];
}

REAL Kinematics::inv_eps(REAL eps,signed int sign) {
  REAL eps2_1 = sqr(eps)-1;
  return 0.5*(Mt2-mQt2+eps2_1+sign*eps*sqrt(sqr(Mt2+mQt2-eps2_1)-4.0*Mt2*mQt2))/eps2_1;
}

void Kinematics::kin1(REAL& z,REAL& l)
{
  //  cerr << "pq = " << pq << ", sigma*z = " << sigma_l << endl;
  REAL zmin = 0;
  REAL zmax = 1;
  if ( eps < eps_max_0 ) {
    zmin = inv_eps(eps,-1);
    if ( eps < eps_max_1 ) {
      zmax = inv_eps(eps,+1);
    }
  }
  //  cerr << zmin << "  " << zmax << endl;  
  FragmentationFunction f(zmin,zmax,Mt2*pq2);
  z = f.getValue();
  l =  get_eps(z);
}

void Kinematics::kin2(REAL& z,REAL& l)
{
  REAL e = rand_gen::Random(eps_min,min(eps_max_0,eps));
  //  if ( e > eps_max_1 )
    z = inv_eps(e,-1);
    //  else
    //    z = inv_eps(e,+1);
  l = e*pq-p0[0];
  if ( get_eps(z) != l ) {
    cerr << "get_eps != eps :  " << get_eps(z) << "  " << l << endl;
    throw "Not possible...";
  }
}

void Kinematics::kin3(REAL& z,REAL& l)
{
  REAL zmin = 0;
  REAL zmax = 1;
  if ( eps < eps_max_0 ) {
    zmin = inv_eps(eps,-1);
    if ( eps < eps_max_1 ) {
      zmax = inv_eps(eps,+1);
    }
  }
  //  cerr << zmin << "  " << zmax << endl;  
  Banerjee f(zmin,zmax);
  z = f.getValue();
  l =  get_eps(z);
}

