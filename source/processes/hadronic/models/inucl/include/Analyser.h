#define WITH_NUCLEI

#ifndef ANALYSER_H
#define ANALYSER_H

#include "CollisionOutput.h"
#include "InuclParticle.h"
#include "InuclElementaryParticle.h"
#include "InuclNuclei.h"
#include "NuclWatcher.h"
#include "ExitonConfiguration.h"

#include "vector"

class Analyser {

public:

Analyser() {

 eventNumber = 0.;
 averageMultiplicity = 0.;
 averageNucleonKinEnergy = 0.;
 averageProtonKinEnergy = 0.;
 averageNeutronKinEnergy = 0.;
 averagePionKinEnergy = 0.;
 averageProtonNumber = 0.;
 averageNeutronNumber = 0.;
 averagePionNumber = 0.;
 averageExitationEnergy = 0.;
 averageNucleiFragments = 0.;
 averagePionPl = 0.;
 averagePionMin = 0.;
 averagePion0 = 0.;
 averageA = 0.;
 averageZ = 0.;
 
 withNuclei = false;
 fissy_prob = 0.;
 
};

void setInelCsec(double csec, bool withn);
 
void setWatchers(const vector<NuclWatcher>& watchers) {
  ana_watchers = watchers;
  cout << " watchers set " << watchers.size() << endl;
};

void try_watchers(double a, double z, bool if_nucl) {
  for(int iw = 0; iw < ana_watchers.size(); iw++) { 
    if(if_nucl) {
      if(ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a,z); 
    }
     else {
      if(!ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a,z); 
    }; 
  };
};

void analyse(const CollisionOutput& output) {
  if(withNuclei) {
    vector<InuclNuclei> nucleus = output.getNucleiFragments();
    if(nucleus.size() >= 0) {
      int nbig = 0;
      averageNucleiFragments += nucleus.size();
      for(int in = 0; in < nucleus.size(); in++) {
        averageExitationEnergy += nucleus[in].getExitationEnergy();
        double a = nucleus[in].getA();
        double z = nucleus[in].getZ();
	if(in == 0) { averageA += a; averageZ += z; };
        if(a > 10.) nbig++;
	try_watchers(a,z,true);
      };
      if(nbig > 1) fissy_prob += 1.;
      eventNumber += 1.;
      vector<InuclElementaryParticle> particles = output.getOutgoingParticles();
      averageMultiplicity += particles.size();
      for(int i = 0; i < particles.size(); i++) {
        double ap;
	double zp;
        if(particles[i].nucleon()) {
	  averageNucleonKinEnergy += particles[i].getKineticEnergy();
          if(particles[i].type() == 1) {
            zp = 1.;
	    ap = 1.;
	    averageProtonNumber += 1.;
	    averageProtonKinEnergy += particles[i].getKineticEnergy();
          }
           else {
	    ap = 1.;
	    zp = 0.;
            averageNeutronNumber += 1.;
	    averageNeutronKinEnergy += particles[i].getKineticEnergy();
          };  
        }
         else if(particles[i].pion()) {
          averagePionKinEnergy += particles[i].getKineticEnergy();
          averagePionNumber += 1.;
          ap = 0.;
	  if(particles[i].type() == 3) {
            zp = 1.;
            averagePionPl += 1.;
	  }
	   else if(particles[i].type() == 5) {  
            zp = -1.;
            averagePionMin += 1.;
	  }
	   else if(particles[i].type() == 7) { 
            zp = 0.;
            averagePion0 += 1.;
	  };
	};
	try_watchers(ap,zp,false);
      };
    };
  }
   else {
      eventNumber += 1.;
      vector<InuclElementaryParticle> particles = output.getOutgoingParticles();
      averageMultiplicity += particles.size();
      for(int i = 0; i < particles.size(); i++) {
        if(particles[i].nucleon()) {
          averageNucleonKinEnergy += particles[i].getKineticEnergy();
          if(particles[i].type() == 1) {
            averageProtonNumber += 1.;
	    averageProtonKinEnergy += particles[i].getKineticEnergy();
          }
           else {
            averageNeutronNumber += 1.;
	    averageNeutronKinEnergy += particles[i].getKineticEnergy();
          };  
        }
         else if(particles[i].pion()) {
          averagePionKinEnergy += particles[i].getKineticEnergy();
          averagePionNumber += 1.;
        };
      };
  }; 
};

void printResults() {
  cout << " Number of events " << int(eventNumber + 0.1) << endl
    << " average multiplicity " << averageMultiplicity/eventNumber << endl
    << " average proton number " << averageProtonNumber/eventNumber << endl
    << " average neutron number " << averageNeutronNumber/eventNumber << endl
    << " average nucleon Ekin " << averageNucleonKinEnergy/
      (averageProtonNumber + averageNeutronNumber) << endl
    << " average proton Ekin " << averageProtonKinEnergy/(averageProtonNumber +
         1.e-10) << endl
    << " average neutron Ekin " << averageNeutronKinEnergy/(averageNeutronNumber +
         1.e-10) << endl
    << " average pion number " << averagePionNumber/eventNumber << endl
    << " average pion Ekin " << averagePionKinEnergy/(averagePionNumber +
                   1.e-10) << endl
    << " average pi+ " << averagePionPl/eventNumber << endl
    << " average pi- " << averagePionMin/eventNumber << endl
    << " average pi0 " << averagePion0/eventNumber << endl;
     		   
  if(withNuclei) {
    cout
    << " average A " << averageA/eventNumber << endl 		   
    << " average Z " << averageZ/eventNumber << endl 		   
    << " average Exitation Energy " << 
         averageExitationEnergy/averageNucleiFragments << endl
    << " average num of fragments " << averageNucleiFragments/eventNumber << endl;
    cout << " fission prob. " << fissy_prob/eventNumber << " c.sec " <<
      inel_csec*fissy_prob/eventNumber << endl;
    handleWatcherStatistics();
  };
};

void printResultsSimple() {
  cout << " Number of events " << int(eventNumber + 0.1) << endl
    << " average multiplicity " << averageMultiplicity/eventNumber << endl
    << " average proton number " << averageProtonNumber/eventNumber << endl
    << " average neutron number " << averageNeutronNumber/eventNumber << endl
    << " average nucleon Ekin " << averageNucleonKinEnergy/
      (averageProtonNumber + averageNeutronNumber) << endl
    << " average proton Ekin " << averageProtonKinEnergy/(averageProtonNumber +
         1.e-10) << endl
    << " average neutron Ekin " << averageNeutronKinEnergy/(averageNeutronNumber +
         1.e-10) << endl
    << " average pion number " << averagePionNumber/eventNumber << endl
    << " average pion Ekin " << averagePionKinEnergy/(averagePionNumber +
                   1.e-10) << endl;
  if(withNuclei) {
    cout		   
    << " average Exitation Energy " << 
         averageExitationEnergy/averageNucleiFragments << endl
    << " average num of fragments " << averageNucleiFragments/eventNumber << endl;
    cout << " fission prob. " << fissy_prob/eventNumber << " c.sec " <<
      inel_csec*fissy_prob/eventNumber << endl;
  };
};

void handleWatcherStatistics();

private: 

double eventNumber;
double averageMultiplicity;
double averageProtonNumber;
double averageNeutronNumber;
double averagePionNumber;
double averageNucleonKinEnergy;
double averageProtonKinEnergy;
double averageNeutronKinEnergy;
double averagePionKinEnergy;
double averageExitationEnergy;
double averageNucleiFragments;
double fissy_prob;
double averagePionPl;
double averagePionMin;
double averagePion0;
double averageA;
double averageZ;

vector<NuclWatcher> ana_watchers;
double inel_csec;
bool withNuclei;

};        

#endif // ANALYSER_H 
