#include "CascadParticle.h"

double CascadParticle::getPathToTheNextZone(double rz_in, double rz_out) {

  double path = -1.;
  double rp = 0.;
  double rr = 0.;
  double pp = 0.;
  vector<double> mom = theParticle.getMomentum();
  for(int i = 1; i < 4; i++) {
    rp += mom[i]*position[i-1];
    rr += position[i-1]*position[i-1];
    pp += mom[i]*mom[i];
  };
  double ra = rr - rp*rp/pp;
  pp = sqrt(pp);
  double ds;
  double d2;
  if(current_zone == 0 || rp > 0.) {
    d2 = rz_out*rz_out - ra;
    ds = 1.;
    movingIn = false;  
  } 
   else { 
    d2 = rz_in*rz_in - ra;
    if(d2 > 0.) {
      ds = -1.;
      movingIn = true;  
    }
     else {
      d2 = rz_out*rz_out - ra;
      ds = 1.;
      movingIn = false;
    };
   };
   path = ds*sqrt(d2) - rp/pp;

   return path;    


}

void CascadParticle::propagateAlongThePath(double path) {
  vector<double> mom = theParticle.getMomentum();
  double pmod = theParticle.getMomModule();
  for(int i = 0; i < 3; i++) position[i] += mom[i+1]*path/pmod;
}

