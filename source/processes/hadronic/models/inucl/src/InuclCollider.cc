//#define DEBUG

#include "InuclCollider.h"
#include "InuclElementaryParticle.h"
#include "LorentzConvertor.h"
#include "ParticleLargerEkin.h"
#include "algorithm"

typedef vector<InuclElementaryParticle>::iterator particleIterator;
typedef vector<InuclNuclei>::iterator nucleiIterator;
	 
CollisionOutput InuclCollider::collide(InuclParticle* bullet,
                     InuclParticle* target) {

  const int itry_max = 1000;
  		     
  CollisionOutput globalOutput;
  InuclElementaryParticle* particle1 =
               dynamic_cast<InuclElementaryParticle*>(bullet);
  InuclElementaryParticle* particle2 =
               dynamic_cast<InuclElementaryParticle*>(target);
  
  if(particle1 && particle2) { // particle + particle 
#ifdef DEBUG
    particle1->printParticle();
    particle2->printParticle();
#endif  
    globalOutput = theElementaryParticleCollider->collide(bullet,target);

  }
   else { // needs to call all machinery    	

     LorentzConvertor convertToTargetRestFrame;

     InteractionCase interCase = bulletTargetSetter(bullet,target);
     int intcase = interCase.getInterCase();
     
     if(intcase > 0) { // ok
       InuclNuclei* ntarget =
                   dynamic_cast<InuclNuclei*>(interCase.getTarget());
       convertToTargetRestFrame.setTarget(ntarget->getMomentum(),
	                                           ntarget->getMass());
       int btype;
       double ab;
       double zb;
       double at = ntarget->getA();
       double zt = ntarget->getZ();
       
       if(intcase == 1) { // particle with nuclei
         InuclElementaryParticle* pbullet = 
	   dynamic_cast<InuclElementaryParticle*>(interCase.getBullet());
         
	 if(pbullet->photon()) {
	   cout << " InuclCollider -> can not collide with photon " << endl;
	   globalOutput.trivialise(bullet,target);
	   return globalOutput;
	 }
	  else {
           convertToTargetRestFrame.setBullet(pbullet->getMomentum(),
	                           pbullet->getMass());   
	   btype = pbullet->type();
	 }; 
       }
        else { // nuclei with nuclei
         InuclNuclei* nbullet = 
	       dynamic_cast<InuclNuclei*>(interCase.getBullet());
         convertToTargetRestFrame.setBullet(nbullet->getMomentum(),
	                           nbullet->getMass());   
         ab = nbullet->getA();
	 zb = nbullet->getZ();
       };
       	
       double ekin = convertToTargetRestFrame.getKinEnergyInTheTRS();
#ifdef DEBUG
       cout << " ekin in trs " << ekin << endl;
#endif
       if(inelasticInteractionPossible(bullet,target,ekin)) {
         convertToTargetRestFrame.toTheTargetRestFrame();
#ifdef DEBUG
	 cout << " degenerated? " << convertToTargetRestFrame.trivial() << endl;
#endif
	 vector<double> bmom(4,0.);
	 bmom[3] = convertToTargetRestFrame.getTRSMomentum();
	 InuclNuclei ntarget(at,zt);
	 vector<double> tmom(4,0.);
	 ntarget.setMomentum(tmom);
	 ntarget.setEnergy();
         theIntraNucleiCascader->setInteractionCase(intcase);
	 
         bool bad = true;
	 int itry = 0;
	 
	 while (bad && itry < itry_max) {
           CollisionOutput TRFoutput;
           CollisionOutput output;
	   itry++;
	   if(intcase == 1) {
	     InuclElementaryParticle pbullet(bmom,btype);
	     output = theIntraNucleiCascader->collide(&pbullet,&ntarget);
	   }
	    else {
	     InuclNuclei nbullet(ab,zb);
	     nbullet.setMomentum(bmom);
	     nbullet.setEnergy();
	     output = theIntraNucleiCascader->collide(&nbullet,&ntarget);
	   };   

#ifdef DEBUG
           cout << " After Cascade " << endl;
           output.printCollisionOutput();
           cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
#endif
//             the rest, if any
           TRFoutput.addOutgoingParticles(output.getOutgoingParticles());
	   if(output.numberOfNucleiFragments() == 1) { // there is smth. after
               
	     InuclNuclei cascad_rec_nuclei = output.getNucleiFragments()[0];

             if(explosion(&cascad_rec_nuclei)) {
               cout << " big bang after cascade " << endl;
	       output = theBigBanger->collide(0,&cascad_rec_nuclei);
               TRFoutput.addOutgoingParticles(output.getOutgoingParticles());
	     }
	      else {
               output = theNonEquilibriumEvaporator->collide(0,&cascad_rec_nuclei);
#ifdef DEBUG
              cout << " After NonEquilibriumEvaporator " << endl;
              output.printCollisionOutput();
              cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
#endif
               TRFoutput.addOutgoingParticles(output.getOutgoingParticles());

               InuclNuclei exiton_rec_nuclei = output.getNucleiFragments()[0];
               output = theEquilibriumEvaporator->collide(0,&exiton_rec_nuclei);
#ifdef DEBUG
               cout << " After EquilibriumEvaporator " << endl;
               output.printCollisionOutput();
               cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
#endif
               TRFoutput.addOutgoingParticles(output.getOutgoingParticles());  
               TRFoutput.addTargetFragments(output.getNucleiFragments());         
	     };
	   };
//             convert to the LAB       
           bool withReflection = convertToTargetRestFrame.reflectionNeeded();       
           vector<InuclElementaryParticle> particles = 
	             TRFoutput.getOutgoingParticles();
	   if(!particles.empty()) { 
	     particleIterator ipart;
	     for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	       vector<double> mom = ipart->getMomentum();
	       if(withReflection) mom[3] = -mom[3];
	       mom = convertToTargetRestFrame.rotate(mom);
	       ipart->setMomentum(mom); 
	       mom = convertToTargetRestFrame.backToTheLab(ipart->getMomentum());
	       ipart->setMomentum(mom); 
	     };
	     sort(particles.begin(), particles.end(), ParticleLargerEkin());
           };
           
	   vector<InuclNuclei> nucleus = 
	             TRFoutput.getNucleiFragments();
	   if(!nucleus.empty()) { 
	     nucleiIterator inuc;
	     for(inuc = nucleus.begin(); inuc != nucleus.end(); inuc++) {
	       vector<double> mom = inuc->getMomentum(); 
	       if(withReflection) mom[3] = -mom[3];
	       mom = convertToTargetRestFrame.rotate(mom);
	       inuc->setMomentum(mom);
	       inuc->setEnergy(); 
	       mom = convertToTargetRestFrame.backToTheLab(inuc->getMomentum());
	       inuc->setMomentum(mom);
	       inuc->setEnergy(); 
	     };
           };
           globalOutput.addOutgoingParticles(particles);
	   globalOutput.addTargetFragments(nucleus);
	   globalOutput.setOnShell(bullet,target);
	   if(globalOutput.acceptable()) {
	     return globalOutput;
           }
	    else {
	     globalOutput.reset();
	   }; 
	 };
#ifdef DEBUG
	 cout << " InuclCollider -> can not generate acceptable inter. after " 
	      << itry_max << " attempts " << endl;
#endif
	 globalOutput.trivialise(bullet,target);
	 return globalOutput;        
       }
        else {
#ifdef DEBUG
	 cout << " InuclCollider -> inelastic interaction is impossible " <<
	   endl << " due to the coulomb barirer " << endl;
#endif
	 globalOutput.trivialise(bullet,target);
	 return globalOutput;
       };	
     }
      else {
       cout << " InuclCollider -> inter case " << intcase << endl;
     };       
  };

  
return globalOutput;

}
		     
bool InuclCollider::inelasticInteractionPossible(InuclParticle* bullet,
          InuclParticle* target, double ekin) const {

const double coeff = 0.001*1.2;
const double one_third = 1./3.;

bool possible = true;
double at;
double zt;
double ab;
double zb;

if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {
  at = nuclei_target->getA();
  zt = nuclei_target->getZ(); 
  if(InuclNuclei* nuclei_bullet = dynamic_cast<InuclNuclei*>(bullet)) {
    ab = nuclei_bullet->getA();
    zb = nuclei_bullet->getZ();     
  }
   else {
    InuclElementaryParticle* particle =
               dynamic_cast<InuclElementaryParticle*>(bullet);
    ab = 1;
    zb = particle->getCharge();
  }; 
}
 else {
  if(InuclNuclei* nuclei_bullet = dynamic_cast<InuclNuclei*>(bullet)) {
    ab = nuclei_bullet->getA();
    zb = nuclei_bullet->getZ();     
    InuclElementaryParticle* particle =
               dynamic_cast<InuclElementaryParticle*>(target);
    at = 1;
    zt = particle->getCharge();    
  }
   else {
    return possible;
  };  
}; 

double VCOL = coeff*zt*zb/(pow(at,one_third) + pow(ab,one_third));
possible = VCOL < ekin;
return possible;

}
	
InteractionCase InuclCollider::bulletTargetSetter(InuclParticle* bullet,
       InuclParticle* target) const {

InteractionCase interCase;

if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {     
  if(InuclNuclei* nuclei_bullet = dynamic_cast<InuclNuclei*>(bullet)) { // A + A         
    interCase.setInterCase(2);
    if(nuclei_target->getA() >= nuclei_bullet->getA()) {
      interCase.setBulletTarget(bullet,target);
    }
     else {
      interCase.setBulletTarget(target,bullet);
    }; 
  }
   else {
    interCase.setInterCase(1);
    interCase.setBulletTarget(bullet,target);
  }; 
}
 else {
  if(InuclNuclei* nuclei_bullet = dynamic_cast<InuclNuclei*>(bullet)) { 
    if(InuclElementaryParticle* part = 
             dynamic_cast<InuclElementaryParticle*>(target)) {
      interCase.setInterCase(1);
      interCase.setBulletTarget(target,bullet);
    };
  };
};
return interCase;

}       

bool InuclCollider::explosion(InuclNuclei* target) const {

const double a_cut = 20.;
const double be_cut = 3.;

  double a = target->getA();
  double z = target->getZ();
  double eexs = target->getExitationEnergy();
  bool explo = true;
  if(a > a_cut) {
    explo = false;
  }
   else {
    if(eexs < be_cut*bindingEnergy(a,z)) explo = false;
  };   
  return explo;

}
 
