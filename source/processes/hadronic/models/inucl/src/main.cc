//#define CHECK_MOMC

#include "Collider.h"
#include "InuclCollider.h"
#include "IntraNucleiCascader.h"
#include "NonEquilibriumEvaporator.h"
#include "EquilibriumEvaporator.h"
#include "Fissioner.h"
#include "BigBanger.h"
#include "ElementaryParticleCollider.h"
#include "InuclParticle.h"
#include "InuclElementaryParticle.h"
#include "InuclNuclei.h"
#include "CollisionOutput.h"
#include "Analyser.h"
#include "WatcherGun.h"

#include "vector"

int main() {
// test program to run inucl stand alone

const int to_report = 100;

int nrain = 100000; // number of interactions to be generated

// auxiliarly stuff for ugly analysis
Analyser* analyser = new Analyser;
WatcherGun* gun = new WatcherGun;
gun->setWatchers();
analyser->setWatchers(gun->getWatchers());
analyser->setInelCsec(1760.,true);

//  colliders initialisation
ElementaryParticleCollider* colep = new ElementaryParticleCollider;
IntraNucleiCascader* incas = new IntraNucleiCascader;
NonEquilibriumEvaporator* noneq = new NonEquilibriumEvaporator;
EquilibriumEvaporator* eqil = new EquilibriumEvaporator;
Fissioner* fiss = new Fissioner;
BigBanger* bigb = new BigBanger;
InuclCollider* inucl = new InuclCollider(colep,incas,noneq,eqil,fiss,bigb);

// Bullet / Target initialisation

// Bullet could be nucleon or pion or nuclei

// 0.8 GeV proton with momentum along Z axis
InuclParticle* bull = new InuclElementaryParticle(0.8,1); 

/*
// neutron  with the momentum defined by the vector momb(4)
vector<double> momb(4);
momb[1] = 0.; momb[2] = 0.; momb[3] = 2.;
InuclParticle* bull = new InuclElementaryParticle(mom,2);
*/

/*
// He4 nuclei with momentum momb
InuclParticle* bull = new InuclNuclei(momb,4.,2.);
((InuclNuclei*)bull)->setEnergy(); // for nuclei mom[4] has to be set specially !!!
*/
 
// Target could be nucleon or nuclei (more precise, in case when
// particle-particle interaction will be generated, at least one particle
// has to be nucleon)

// Au197 target at rest
InuclParticle* targ = new InuclNuclei(0.,197.,79.);

/*
// Neutron with momentum momta
vector<double> momta(4);
momta[1] = -0.3; momta[2] = 0.2; momta[3] = -0.2; 
InuclParticle* targ = new InuclElementaryParticle(momta,2);
*/

/*
// C12 nuclei with momentum momta
InuclParticle* targ = new InuclNuclei(momta,12.,6.);
((InuclNuclei*)targ)->setEnergy();
*/

#ifdef CHECK_MOMC
vector<double> total_mom_in = bull->getMomentum();
vector<double> momt = targ->getMomentum();
for(int i = 0; i < 4; i++) total_mom_in[i] += momt[i];
vector<double> total_mom_out;
bull->printParticle();
targ->printParticle();
cout << " tot in mom: px " << total_mom_in[1] << " py " << total_mom_in[2]
 << " pz " << total_mom_in[3] << " e " << total_mom_in[0] << endl;
#endif

for(int i = 0; i < nrain; i++) {
  if((i+1) % to_report == 0)  
    cout << " ---- EVENT ------ " << i+1 << endl;
    CollisionOutput output = inucl->collide(bull,targ); // standard method
//  auxiliarly method to get more information about the different stages 
//  of interaction  
//  CollisionOutput output = inucl->collideAndTest(bull,targ,analyser);
//
//  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
//  cout << " Global output " << endl;
//  output.printCollisionOutput();
//
//
#ifdef CHECK_MOMC
  total_mom_out = output.getTotalOutputMomentum();
  cout << " 4 - momentum conservation check " << endl
    << " dE " << total_mom_out[0] - total_mom_in[0] << 
     " dPx " << total_mom_out[1] - total_mom_in[1] <<
     " dPy " << total_mom_out[2] - total_mom_in[2] <<
     " dPz " << total_mom_out[3] - total_mom_in[3] << endl;
#endif
  analyser->analyse(output);
};
analyser->printResults();

return 0;
        
}   
