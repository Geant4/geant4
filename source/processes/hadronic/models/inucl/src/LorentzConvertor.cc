//#define DEBUG

#include "LorentzConvertor.h"

void LorentzConvertor::toTheCenterOfMass() {
   
const double small = 1.e-10;

v2 = 0.;
double pv = 0.;
double e_sum = target_mom[0]  + bullet_mom[0];
velocity.resize(4);
for(int i = 1; i < 4; i++) {
  velocity[i] = (target_mom[i]  + bullet_mom[i])/e_sum;
  v2 += velocity[i]*velocity[i];
  pv += target_mom[i]*velocity[i];
};
   
gamma = 1./sqrt(fabs(1. - v2));
ecm_tot = e_sum/gamma;

double pa = 0.;
double pb = 0.;
scm_momentum.resize(4);
double xx = pv*(gamma-1.)/v2 - target_mom[0]*gamma;
for(int i = 1; i < 4; i++) {
  scm_momentum[i] = -target_mom[i] -velocity[i]*xx;
#ifdef DEBUG
  cout << " i " << i << " pscm(i) " << scm_momentum[i] << endl;
#endif
  pa += scm_momentum[i]*scm_momentum[i];
  pb += scm_momentum[i]*velocity[i];
};
ga = v2 - pb*pb/pa;
if(ga < small) {
  ga = small;
  degenerated = true;
#ifdef DEBUG
  cout << " degenerated case " << endl; 
#endif
}
 else {
  ga = sqrt(ga);
}; 
#ifdef DEBUG
cout << " ga " << ga << " v2 " << v2 << " pb " << pb << " pb*pb/pa " << pb*pb/pa
<< " pv " << pv << endl;
#endif
pscm = sqrt(pa);
gb = pb/pscm;
gbpp = gb/pscm;
gapp = ga*pscm;
   
}

vector<double> LorentzConvertor::rotate(const vector<double> mom) const {
vector<double> mom_rot(4);

#ifdef DEBUG
cout << " ga " << ga << " gbpp " << gbpp << " gapp " << gapp << endl;  
cout << " gegenerated " << degenerated << endl;
cout << " before rotation: px " << mom[1] << " py " << mom[2] <<
  " pz " << mom[3] << endl;
#endif

if(degenerated) {
  mom_rot = mom; 
}
 else {
  mom_rot[1] = mom[1]*(velocity[1] - gbpp*scm_momentum[1])/ga + 
   mom[2]*(scm_momentum[2]*velocity[3] - scm_momentum[3]*velocity[2])/gapp +
   mom[3]*scm_momentum[1]/pscm;
  mom_rot[2] = mom[1]*(velocity[2] - gbpp*scm_momentum[2])/ga + 
   mom[2]*(scm_momentum[3]*velocity[1] - scm_momentum[1]*velocity[3])/gapp +
   mom[3]*scm_momentum[2]/pscm;
  mom_rot[3] = mom[1]*(velocity[3] - gbpp*scm_momentum[3])/ga + 
   mom[2]*(scm_momentum[1]*velocity[2] - scm_momentum[2]*velocity[1])/gapp +
   mom[3]*scm_momentum[3]/pscm;
};

#ifdef DEBUG
cout << " after rotation: px " << mom_rot[1] << " py " << mom_rot[2] <<
  " pz " << mom_rot[3] << endl;
#endif

return mom_rot;

}

vector<double> LorentzConvertor::rotate(const vector<double> mom1, const vector<double> mom) const {

const double small = 1.e-10;

vector<double> mom_rot(4);

double pp = 0.;
double pv = 0.;
for(int i = 0; i < 4; i++) {
  pp += mom1[i]*mom1[i];
  pv += mom1[i]*velocity[i];
};
double ga1 = v2 - pv*pv/pp;
if(ga1 < small) {
  mom_rot = mom;
}
 else {  
  ga1 = sqrt(ga1);
  double gb1 = pv/pp;
  pp = sqrt(pp);
  double ga1pp = ga1*pp;

  mom_rot[1] = mom[1]*(velocity[1] - gb1*mom1[1])/ga1 + 
   mom[2]*(mom1[2]*velocity[3] - mom1[3]*velocity[2])/ga1pp +
   mom[3]*mom1[1]/pp;
  mom_rot[2] = mom[1]*(velocity[2] - gb1*mom1[2])/ga1 + 
   mom[2]*(mom1[3]*velocity[1] - mom1[1]*velocity[3])/ga1pp +
   mom[3]*mom1[2]/pp;
  mom_rot[3] = mom[1]*(velocity[3] - gb1*mom1[3])/ga1 + 
   mom[2]*(mom1[1]*velocity[2] - mom1[2]*velocity[1])/ga1pp +
   mom[3]*mom1[3]/pp;
};
return mom_rot;

}

void LorentzConvertor::toTheTargetRestFrame() {
   
const double small = 1.e-10;

gamma = target_mom[0]/target_mass;
v2 = 0.;
double pv = 0.;
double e_sum = target_mom[0]  + bullet_mom[0];
velocity.resize(4);
for(int i = 1; i < 4; i++) {
  velocity[i] = target_mom[i]/target_mom[0];
  v2 += velocity[i]*velocity[i];
  pv += bullet_mom[i]*velocity[i];
};
double pa = 0.;
double pb = 0.;
scm_momentum.resize(4);
double xx = 0.;
if(v2 > small) xx = 
               pv*(gamma - 1.)/v2 - bullet_mom[0]*gamma;
for(int i = 1; i < 4; i++) {
    scm_momentum[i] = bullet_mom[i]  + velocity[i]*xx;
#ifdef DEBUG
    cout << " rf: i " << i << " pscm(i) " << scm_momentum[i] << endl;
#endif
    pa += scm_momentum[i]*scm_momentum[i];
    pb += scm_momentum[i]*velocity[i];
};

ga = v2 - pb*pb/pa;
if(ga < small) {
    ga = small;
    degenerated = true;
}
 else {
    ga = sqrt(ga);
};  
pscm = sqrt(pa);
plab = pscm;
gb = pb/pscm;
gbpp = gb/pscm;
gapp = ga*pscm;   

}

vector<double> LorentzConvertor::backToTheLab(const vector<double>& mom) const {

const double small = 1.e-10;

#ifdef DEBUG
cout << " at rest: px " << mom[1] << " py " << mom[2] << " pz " << mom[3] << 
  " e " << mom[0] << endl;
cout << " v2 " << v2 << endl;   
#endif

vector<double> mom1(4);
if(v2 < small) {
  mom1 = mom;
}
 else { 
  double pv = 0.;
  for(int i = 1; i < 4; i++) pv += mom[i]*velocity[i];
  double xx = pv*(gamma - 1.)/v2 + mom[0]*gamma;
  for(int i = 1; i < 4; i++) mom1[i] = mom[i] + velocity[i]*xx;
};
#ifdef DEBUG
cout << " at lab: px " << mom1[1] << " py " << mom1[2] << " pz " << mom1[3] << endl;
#endif

return mom1;

}

bool LorentzConvertor::reflectionNeeded() const {
const double small = 1.e-10;

  if(v2 < small) {
    return false;
  }
   else {   
    if(degenerated) return scm_momentum[3] < 0.;
  };
}
