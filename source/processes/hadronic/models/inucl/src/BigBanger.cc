//#define DEBUG

#include "BigBanger.h"
#include "InuclNuclei.h"
#include "ParticleLargerEkin.h"
#include "LorentzConvertor.h"
#include "algorithm"

typedef vector<InuclElementaryParticle>::iterator particleIterator;

CollisionOutput BigBanger::collide(InuclParticle* bullet,
                     InuclParticle* target) {
// primitive explosion model A -> nucleons to prevent too exotic evaporation

const double small_ekin = 1.e-6;

CollisionOutput output;

if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {
  
  double A = nuclei_target->getA();
  double Z = nuclei_target->getZ();
  vector<double> PEX = nuclei_target->getMomentum();
  double EEXS = nuclei_target->getExitationEnergy();

  InuclElementaryParticle dummy(small_ekin,1);
  LorentzConvertor toTheNucleiSystemRestFrame;
  toTheNucleiSystemRestFrame.setBullet(dummy.getMomentum(),
                                               dummy.getMass());
  toTheNucleiSystemRestFrame.setTarget(PEX,nuclei_target->getMass());
  toTheNucleiSystemRestFrame.toTheTargetRestFrame();

  double etot = 0.001*(EEXS - bindingEnergy(A,Z));
#ifdef DEBUG
  cout << " BigBanger: target " << endl;
  nuclei_target->printParticle(); 
  cout << " BigBanger: a " << A << " z " << Z << " eexs " << EEXS << " etot " <<
    etot << " nm " << nuclei_target->getMass() << endl;
#endif
  
  vector<InuclElementaryParticle> particles = 	    
           generateBangInSCM(etot,A,Z,dummy.getParticleMass(1),
	                                    dummy.getParticleMass(2));
#ifdef DEBUG
        cout << " particles " << particles.size() << endl;
	for(int i = 0; i < particles.size(); i++) 
	  particles[i].printParticle();
#endif
  if(!particles.empty()) { // convert back to Lab
#ifdef DEBUG
	  vector<double> totscm(4,0.);
	  vector<double> totlab(4,0.);
#endif
     particleIterator ipart;
     for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
#ifdef DEBUG
	    vector<double> mom_scm = ipart->getMomentum();
	    for(int i = 0; i < 4; i++) totscm[i] += mom_scm[i];
#endif
       vector<double> mom = 
	      toTheNucleiSystemRestFrame.backToTheLab(ipart->getMomentum());
       ipart->setMomentum(mom); 
#ifdef DEBUG
       mom = ipart->getMomentum();
       for(int i = 0; i < 4; i++) totlab[i] += mom[i];
#endif
     };
     sort(particles.begin(), particles.end(), ParticleLargerEkin());
#ifdef DEBUG
	  cout << " In SCM: total outgoing momentum " << endl 
	   << " E " << totscm[0] << " px " << totscm[1]
	    << " py " << totscm[2] << " pz " << totscm[3] << endl; 
	  cout << " In Lab: mom cons " << endl 
	   << " E " << PEX[0] + 0.001*EEXS - totlab[0] 
	   << " px " << PEX[1] - totlab[1]
	   << " py " << PEX[2] - totlab[2] 
	   << " pz " << PEX[3] - totlab[3] << endl; 
#endif
  };
	
  output.addOutgoingParticles(particles);

}
 else {
  cout << " BigBanger -> try to bang not nuclei " << endl;
}; 

return output;

}		     

vector<InuclElementaryParticle>  	    
    BigBanger::generateBangInSCM(double etot, double a, double z, double mp,
	                      double mn) const {
  const double ang_cut = 0.9999;
  const int itry_max = 1000;
  
  int ia = int(a + 0.1);
  int iz = int(z + 0.1);
#ifdef DEBUG
  cout << " ia " << ia << " iz " << iz << endl;
#endif
  vector<InuclElementaryParticle> particles;
  
  if(ia == 1) {
//    abnormal situation
    double m = iz > 0 ? mp : mn;
    double pmod = sqrt((etot + 2.*m)*etot);
    vector<double> mom(4);
    pair<double,double> COS_SIN = randomCOS_SIN();
    double FI = randomPHI();
    double Pt = pmod*COS_SIN.second;
    mom[1] = Pt*cos(FI);
    mom[2] = Pt*sin(FI);
    mom[3] = Pt*COS_SIN.first;    
    int knd = iz > 0 ? 1 : 2;
    particles.push_back(InuclElementaryParticle(mom,knd));
    return particles;
  };  
     
  vector<double> pmod = generateMomentumModules(etot,a,z,mp,mn);
  bool bad = true;
  int itry = 0;

  while(bad && itry < itry_max) {
  
    itry++;
    vector<vector<double> > scm_momentums;
    vector<double> tot_mom(4);

    if(ia == 2) {
      vector<double> mom(4);
      pair<double,double> COS_SIN = randomCOS_SIN();
      double FI = randomPHI();
      double Pt = pmod[0]*COS_SIN.second;
      mom[1] = Pt*cos(FI);
      mom[2] = Pt*sin(FI);
      mom[3] = Pt*COS_SIN.first;    
      for(int j = 1; j < 4; j++) tot_mom[j] += mom[j];		 
      scm_momentums.push_back(mom);
      vector<double> mom1(4);
      for(int i = 1; i < 4; i++) mom1[i] = - mom[i];
      scm_momentums.push_back(mom1);  
      bad = false;
    }
     else {
       for(int i = 0; i < ia - 2; i++) {
         vector<double> mom(4);
         pair<double,double> COS_SIN = randomCOS_SIN();
         double FI = randomPHI();
         double Pt = pmod[i]*COS_SIN.second;
         mom[1] = Pt*cos(FI);
         mom[2] = Pt*sin(FI);
         mom[3] = Pt*COS_SIN.first;    
         for(int j = 1; j < 4; j++) tot_mom[j] += mom[j];		 
         scm_momentums.push_back(mom);
       };

//                handle last two
       double tot_mod = sqrt(tot_mom[1]*tot_mom[1] + 
	         tot_mom[2]*tot_mom[2] + tot_mom[3]*tot_mom[3]); 
       double ct = -0.5*(tot_mod*tot_mod + pmod[ia-2]*pmod[ia-2] -
                 pmod[ia-1]*pmod[ia-1])/tot_mod/pmod[ia-2];

#ifdef DEBUG
	     cout << " ct last " << ct << endl;
#endif
  
       if(fabs(ct) < ang_cut) {
         vector<double> mom2 = 
		     generateWithFixedTheta(ct,pmod[ia-2]);
//       rotate to the normal system
         vector<double> apr = tot_mom;
         for(int i = 1; i < 4; i++) apr[i] /= tot_mod;
         double a_tr = sqrt(apr[1]*apr[1] + apr[2]*apr[2]);
         vector<double> mom(4);
         mom[1] = mom2[3]*apr[1] + (mom2[1]*apr[2] + mom2[2]*apr[3]*apr[1])/a_tr;
         mom[2] = mom2[3]*apr[2] + (-mom2[1]*apr[1] + mom2[2]*apr[3]*apr[2])/a_tr;      
         mom[3] = mom2[3]*apr[3] - mom2[2]*a_tr;      
         scm_momentums.push_back(mom);
//               and the last one
         vector<double> mom1(4);
         for(int i = 1; i < 4; i++) mom1[i] = - mom[i] - tot_mom[i];
         scm_momentums.push_back(mom1);  
         bad = false;
       };
    };   
    if(!bad) {
      for(int i = 0; i < ia; i++) {
         int knd = i < iz ? 1 : 2;
         particles.push_back(InuclElementaryParticle(
		                      scm_momentums[i],knd));
      };
    };
  };  
#ifdef DEBUG
  if(itry == itry_max) cout << " BigBanger-> can not generate bang " << endl;
#endif
  return particles;
  
}
	   
vector<double> BigBanger::generateMomentumModules(double etot, 
            double a, double z, double mp, double mn) const {

  int ia = int(a + 0.1);
  int iz = int(z + 0.1);
  vector<double> pmod;
  double xtot = 0.;
  double promax = maxProbability(a);
  
  for(int i = 0; i < ia; i++) { 
    double x = generateX(ia,a,promax);
#ifdef DEBUG
    cout << " i " << i << " x " << x << endl;
#endif
    pmod.push_back(x);
    xtot += x;
  };
  for(int i = 0; i < ia; i++) {
    double m = i < iz ? mp : mn;
    pmod[i] = pmod[i]*etot/xtot;
    pmod[i] = sqrt(pmod[i]*(pmod[i] + 2.*m));
#ifdef DEBUG
    cout << " i " << i << " pmod " << pmod[i] << endl;
#endif
  };
  return pmod;
  
}

double BigBanger::xProbability(double x, int ia) const {
  int ihalf = ia/2;
  double ekpr = 0.;
  if(x < 1. || x > 0.) {
    ekpr = x*x;
    if(2*ihalf == ia) { // even A
      ekpr *= sqrt(1. - x)*pow((1. -x),int((3*ia - 6)/2)); 
    }
     else {
      ekpr *= pow((1. -x),int((3*ia - 5)/2));
    };
  }; 
  
  return ekpr;

}

double BigBanger::maxProbability(double a) const {
  return xProbability(1./(a - 1.)/1.5,int(a + 0.1));
}

double BigBanger::generateX(int ia, double a, double promax) const {
  const int itry_max = 1000;
  int itry = 0;
  double x;
  
  while(itry < itry_max) {
    itry++;
    x = inuclRndm();
    if(xProbability(x,ia) >= promax*inuclRndm()) return x;
  };
#ifdef DEBUG
  cout << " BigBanger-> can not generate x " << endl;
#endif
  return maxProbability(a);

}
