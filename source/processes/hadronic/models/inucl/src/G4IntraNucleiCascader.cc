//#define DEBUG
#define RUN

#include "G4IntraNucleiCascader.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzConvertor.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4NucleiModel.hh"
#include "G4CascadParticle.hh"
#include "algorithm"


typedef vector<G4InuclElementaryParticle>::iterator particleIterator;

G4CollisionOutput G4IntraNucleiCascader::collide(G4InuclParticle* bullet,
                                 G4InuclParticle* target) {

const G4int itry_max = 1000;
const G4int reflection_cut = 500;
const G4double eexs_cut = 0.0001;

#ifdef DEBUG
bullet->printParticle();
target->printParticle();
#endif

G4CollisionOutput output;

#ifdef RUN

G4InuclNuclei* tnuclei = dynamic_cast<G4InuclNuclei*>(target);
G4InuclNuclei* bnuclei = dynamic_cast<G4InuclNuclei*>(bullet);
G4InuclElementaryParticle* bparticle = dynamic_cast<G4InuclElementaryParticle*>
  (bullet);
NucleiModel model(tnuclei);
vector<double> momentum_in = bullet->getMomentum();
momentum_in[0] += tnuclei->getMass();
double ekin_in; 
#ifdef DEBUG
model.printModel();
cout << " intitial momentum  E " << momentum_in[0] << " Px " << momentum_in[1] 
   << " Py " << momentum_in[2] << " Pz " << momentum_in[3] << endl;
#endif

int itry = 0;

while(itry < itry_max) {
  itry++;
  model.reset();
  vector<CascadParticle> cascad_particles;
  ExitonConfiguration theExitonConfiguration;
  vector<InuclElementaryParticle> output_particles;
  double afin = tnuclei->getA();
  double zfin = tnuclei->getZ();
   
  if(inter_case == 1) { // particle with nuclei
    ekin_in = bparticle->getKineticEnergy();
    zfin += bparticle->getCharge();
    if(bparticle->nucleon()) afin += 1.;
    cascad_particles.push_back(model.initializeCascad(bparticle));
  }
   else { // nuclei with nuclei
    ekin_in = bnuclei->getKineticEnergy();
    double ab = bnuclei->getA();
    double zb = bnuclei->getZ();
    afin += ab;
    zfin += zb;
    pair<vector<CascadParticle>, vector<InuclElementaryParticle> > 
       all_particles = model.initializeCascad(bnuclei,tnuclei);
    cascad_particles = all_particles.first;
    for(int ip = 0; ip < all_particles.second.size(); ip++) 
      output_particles.push_back(all_particles.second[ip]);
    if(cascad_particles.size() == 0) { // compound nuclei
      int ia = int(ab + 0.5);
      int iz = int(zb + 0.5);
      for(int i = 0; i < ia; i++) {
        int knd = i < iz ? 1 : 2;
	theExitonConfiguration.incrementQP(knd);
      };
      int ihn = int(2.*(ab - zb)*inuclRndm() + 0.5);
      int ihz = int(2.*zb*inuclRndm() + 0.5);
      for(int i = 0; i < ihn; i++) theExitonConfiguration.incrementHoles(2);
      for(int i = 0; i < ihz; i++) theExitonConfiguration.incrementHoles(1);
    };
  }; 

  vector<CascadParticle> new_cascad_particles;
  int iloop = 0;
  while(!cascad_particles.empty() && !model.empty()) {
    iloop++;
#ifdef DEBUG
    cout << " ***** number of cparticles " << cascad_particles.size() << endl;
    cascad_particles.back().print();
#endif

    new_cascad_particles = model.generateParticleFate(cascad_particles.back(),
                                             theElementaryParticleCollider);

#ifdef DEBUG
    cout << " new particles " << new_cascad_particles.size() << endl;
#endif
//  handle the result of a new step
    if(new_cascad_particles.size() == 1) { // last particle goes without interaction
      cascad_particles.pop_back();
      if(model.stillInside(new_cascad_particles[0])) { // particle survives 
#ifdef DEBUG
        cout << " still inside " << endl;
#endif
        if(new_cascad_particles[0].getNumberOfReflections() < reflection_cut &&
	   model.worthToPropagate(new_cascad_particles[0])) { // it's ok
#ifdef DEBUG
          cout << " survives " << endl;
#endif
	  cascad_particles.push_back(new_cascad_particles[0]);
	}
	 else { // it becomes an exiton 
#ifdef DEBUG
          cout << " becomes an exiton " << endl;
#endif
	  theExitonConfiguration.incrementQP(new_cascad_particles[0].getParticle().type());
	};  
      }
       else { // goes out
#ifdef DEBUG
          cout << " **** goes out **** " << endl;
	  new_cascad_particles[0].print();
#endif
        output_particles.push_back(new_cascad_particles[0].getParticle());
      }; 
    }
     else { // interaction 
      cascad_particles.pop_back();
      for(int i = 0; i < new_cascad_particles.size(); i++) 
                           cascad_particles.push_back(new_cascad_particles[i]);
      pair<int,int> holes = model.getTypesOfNucleonsInvolved();
      theExitonConfiguration.incrementHoles(holes.first);
      if(holes.second > 0) theExitonConfiguration.incrementHoles(holes.second);
    };
  };

//  cascad is finished -> check, whether it's o'k
#ifdef DEBUG
      cout << " ***** cascad finished ******* " << endl
       << " output_particles  " << output_particles.size() <<  endl;
#endif
  vector<double> momentum_out(4,0.);
  particleIterator ipart;
  for(ipart = output_particles.begin(); ipart != output_particles.end(); ipart++) {
    vector<double> mom = ipart->getMomentum();
    for(int j = 0; j < 4; j++) momentum_out[j] += mom[j];
    zfin -= ipart->getCharge();
    if(ipart->nucleon()) afin -= 1.;
  };
#ifdef DEBUG
      cout << "  afin " << afin << " zfin " << zfin <<  endl;
#endif
  if(afin > 1.) {
    InuclNuclei outgoing_nuclei(afin,zfin);
    double mass = outgoing_nuclei.getMass();
    momentum_out[0] += mass;        
    for(int j = 0; j < 4; j++) momentum_out[j] = momentum_in[j] - momentum_out[j];
#ifdef DEBUG
      cout << "  Eex + Ekin " << momentum_out[0]  <<  endl;
#endif
    if(momentum_out[0] > 0.) { // Eex + Ekin > 0.
      double pnuc = momentum_out[1]*momentum_out[1] + momentum_out[2]*momentum_out[2] +
             momentum_out[3]*momentum_out[3]; 
      double ekin = sqrt(mass*mass + pnuc) - mass;
      double Eex = 1000.*(momentum_out[0] - ekin);
#ifdef DEBUG
      cout << "  Eex  " << Eex  <<  endl;
#endif
      if(goodCase(afin,zfin,Eex,ekin_in)) { // ok, exitation energy > cut
        sort(output_particles.begin(), output_particles.end(), ParticleLargerEkin());
        output.addOutgoingParticles(output_particles);
        outgoing_nuclei.setMomentum(momentum_out);
        outgoing_nuclei.setEnergy();
        outgoing_nuclei.setExitationEnergy(Eex);
        outgoing_nuclei.setExitonConfiguration(theExitonConfiguration);	                           
        output.addTargetFragment(outgoing_nuclei);
        return output;
      };
    };
  }
   else { // special case, when one has no nuclei after the cascad
    if(afin == 1.) { // recoiling nucleon
      for(int j = 0; j < 4; j++) momentum_out[j] = momentum_in[j] - momentum_out[j];
      InuclElementaryParticle  last_particle;
      if(zfin == 1.) { // recoiling proton
        last_particle.setType(1);
      }
       else { // neutron
        last_particle.setType(2);
      }; 
      last_particle.setMomentum(momentum_out);
      output_particles.push_back(last_particle);
    }; 
    sort(output_particles.begin(), output_particles.end(), ParticleLargerEkin());
    output.addOutgoingParticles(output_particles);
    return output;
  }; 
};

#else
//  special branch to avoid the cascad generation but to get the input for
//  evaporation etc
//
  vector<double> momentum_out(4,0.);
  InuclNuclei outgoing_nuclei(169,69);
  outgoing_nuclei.setMomentum(momentum_out);
  outgoing_nuclei.setEnergy();
  outgoing_nuclei.setExitationEnergy(150.);
  ExitonConfiguration theExitonConfiguration(3.,3.,5.,6.);
  outgoing_nuclei.setExitonConfiguration(theExitonConfiguration);	                           
  output.addTargetFragment(outgoing_nuclei);
  return output;
//
/*
InuclElementaryParticle* bparticle = dynamic_cast<InuclElementaryParticle*>
  (bullet);
InuclNuclei* tnuclei = dynamic_cast<InuclNuclei*>(target);
  output.addOutgoingParticle(*bparticle);
  output.addTargetFragment(*tnuclei);
  return output;
*/
#endif

#ifdef DEBUG
cout << " IntraNucleiCascader-> no inelastic interaction after " << itry_max << " attempts "
     << endl;
#endif
output.trivialise(bullet,target);
return output;

}

bool IntraNucleiCascader::goodCase(double a, double z, double eexs, double ein)
               const {
const double eexs_cut = 0.0001;
const double reason_cut = 7.;
const double ediv_cut = 5.;

bool good = false;

if(eexs > eexs_cut) {
  double eexs_max0z = 1000.*ein/ediv_cut;
  double dm = bindingEnergy(a,z);
  double eexs_max = eexs_max0z > reason_cut*dm ? eexs_max0z : reason_cut*dm;
  if(eexs < eexs_max) good = true;
#ifdef DEBUG
  cout << " eexs " << eexs << " max " << eexs_max << " dm " << dm << endl;
#endif
};

return good; 

}
