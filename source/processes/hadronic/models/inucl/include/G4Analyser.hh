#define WITH_NUCLEI

#ifndef G4ANALYSER_HH
#define G4ANALYSER_HH

#include "G4CollisionOutput.hh"
#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4NuclWatcher.hh"
#include "G4ExitonConfiguration.hh"

#include "vector"

class G4Analyser {

public:

G4Analyser() {

 eventNumber = 0.0;
 averageMultiplicity = 0.0;
 averageNucleonKinEnergy = 0.0;
 averageProtonKinEnergy = 0.0;
 averageNeutronKinEnergy = 0.0;
 averagePionKinEnergy = 0.0;
 averageProtonNumber = 0.0;
 averageNeutronNumber = 0.0;
 averagePionNumber = 0.0;
 averageExitationEnergy = 0.0;
 averageNucleiFragments = 0.0;
 averagePionPl = 0.0;
 averagePionMin = 0.0;
 averagePion0 = 0.0;
 averageA = 0.0;
 averageZ = 0.0;
 
 withNuclei = false;
 fissy_prob = 0.0;
 
};

void setInelCsec(G4double csec, bool withn);
 
void setWatchers(const vector<G4NuclWatcher>& watchers) {
  ana_watchers = watchers;
  G4cout << " watchers set " << watchers.size() << G4endl;
};

void try_watchers(G4double a, G4double z, bool if_nucl) {
  for(G4int iw = 0; iw < ana_watchers.size(); iw++) { 
    if(if_nucl) {
      if(ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a,z); 
    }
     else {
      if(!ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a,z); 
    }; 
  };
};

void analyse(const G4CollisionOutput& output) {
  if(withNuclei) {
    vector<G4InuclNuclei> nucleus = output.getNucleiFragments();
    if(nucleus.size() >= 0) {
      G4int nbig = 0;
      averageNucleiFragments += nucleus.size();
      for(G4int in = 0; in < nucleus.size(); in++) {
        averageExitationEnergy += nucleus[in].getExitationEnergy();
        G4double a = nucleus[in].getA();
        G4double z = nucleus[in].getZ();
	if(in == 0) { averageA += a; averageZ += z; };
        if(a > 10.0) nbig++;
	try_watchers(a, z, true);
      };
      if(nbig > 1) fissy_prob += 1.0;
      eventNumber += 1.0;
      vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();
      averageMultiplicity += particles.size();
      for(G4int i = 0; i < particles.size(); i++) {
        G4double ap;
	G4double zp;
        if(particles[i].nucleon()) {
	  averageNucleonKinEnergy += particles[i].getKineticEnergy();
          if(particles[i].type() == 1) {
            zp = 1.;
	    ap = 1.;
	    averageProtonNumber += 1.0;
	    averageProtonKinEnergy += particles[i].getKineticEnergy();
          }
           else {
	    ap = 1.0;
	    zp = 0.0;
            averageNeutronNumber += 1.0;
	    averageNeutronKinEnergy += particles[i].getKineticEnergy();
          };  
        }
         else if(particles[i].pion()) {
          averagePionKinEnergy += particles[i].getKineticEnergy();
          averagePionNumber += 1.0;
          ap = 0.0;
	  if(particles[i].type() == 3) {
            zp = 1.0;
            averagePionPl += 1.0;
	  }
	   else if(particles[i].type() == 5) {  
            zp = -1.0;
            averagePionMin += 1.0;
	  }
	   else if(particles[i].type() == 7) { 
            zp = 0.0;
            averagePion0 += 1.0;
	  };
	};
	try_watchers(ap, zp, false);
      };
    };
  }
   else {
      eventNumber += 1.0;
      vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();
      averageMultiplicity += particles.size();
      for(G4int i = 0; i < particles.size(); i++) {
        if(particles[i].nucleon()) {
          averageNucleonKinEnergy += particles[i].getKineticEnergy();
          if(particles[i].type() == 1) {
            averageProtonNumber += 1.0;
	    averageProtonKinEnergy += particles[i].getKineticEnergy();
          }
           else {
            averageNeutronNumber += 1.0;
	    averageNeutronKinEnergy += particles[i].getKineticEnergy();
          };  
        }
         else if(particles[i].pion()) {
          averagePionKinEnergy += particles[i].getKineticEnergy();
          averagePionNumber += 1.0;
        };
      };
  }; 
};

void printResults() {
  G4cout << " Number of events " << int(eventNumber + 0.1) << G4endl
    << " average multiplicity " << averageMultiplicity/eventNumber << G4endl
    << " average proton number " << averageProtonNumber/eventNumber << G4endl
    << " average neutron number " << averageNeutronNumber/eventNumber << G4endl
    << " average nucleon Ekin " << averageNucleonKinEnergy/
      (averageProtonNumber + averageNeutronNumber) << G4endl
    << " average proton Ekin " << averageProtonKinEnergy/(averageProtonNumber +
         1.e-10) << G4endl
    << " average neutron Ekin " << averageNeutronKinEnergy/(averageNeutronNumber +
         1.e-10) << G4endl
    << " average pion number " << averagePionNumber/eventNumber << G4endl
    << " average pion Ekin " << averagePionKinEnergy/(averagePionNumber +
                   1.e-10) << G4endl
    << " average pi+ " << averagePionPl/eventNumber << G4endl
    << " average pi- " << averagePionMin/eventNumber << G4endl
    << " average pi0 " << averagePion0/eventNumber << G4endl;
     		   
  if(withNuclei) {
    G4cout
    << " average A " << averageA/eventNumber << G4endl 		   
    << " average Z " << averageZ/eventNumber << G4endl 		   
    << " average Exitation Energy " << 
         averageExitationEnergy/averageNucleiFragments << G4endl
    << " average num of fragments " << averageNucleiFragments/eventNumber << G4endl;
    G4cout << " fission prob. " << fissy_prob/eventNumber << " c.sec " <<
      inel_csec*fissy_prob/eventNumber << G4endl;
    handleWatcherStatistics();
  };
};

void printResultsSimple() {
  G4cout << " Number of events " << int(eventNumber + 0.1) << G4endl
    << " average multiplicity " << averageMultiplicity/eventNumber << G4endl
    << " average proton number " << averageProtonNumber/eventNumber << G4endl
    << " average neutron number " << averageNeutronNumber/eventNumber << G4endl
    << " average nucleon Ekin " << averageNucleonKinEnergy/
      (averageProtonNumber + averageNeutronNumber) << G4endl
    << " average proton Ekin " << averageProtonKinEnergy/(averageProtonNumber +
         1.e-10) << G4endl
    << " average neutron Ekin " << averageNeutronKinEnergy/(averageNeutronNumber +
         1.e-10) << G4endl
    << " average pion number " << averagePionNumber/eventNumber << G4endl
    << " average pion Ekin " << averagePionKinEnergy/(averagePionNumber +
                   1.e-10) << G4endl;
  if(withNuclei) {
    G4cout		   
    << " average Exitation Energy " << 
         averageExitationEnergy/averageNucleiFragments << G4endl
    << " average num of fragments " << averageNucleiFragments/eventNumber << G4endl;
    G4cout << " fission prob. " << fissy_prob/eventNumber << " c.sec " <<
      inel_csec*fissy_prob/eventNumber << G4endl;
  };
};

void handleWatcherStatistics();

private: 

G4double eventNumber;
G4double averageMultiplicity;
G4double averageProtonNumber;
G4double averageNeutronNumber;
G4double averagePionNumber;
G4double averageNucleonKinEnergy;
G4double averageProtonKinEnergy;
G4double averageNeutronKinEnergy;
G4double averagePionKinEnergy;
G4double averageExitationEnergy;
G4double averageNucleiFragments;
G4double fissy_prob;
G4double averagePionPl;
G4double averagePionMin;
G4double averagePion0;
G4double averageA;
G4double averageZ;

vector<G4NuclWatcher> ana_watchers;
G4double inel_csec;
bool withNuclei;

};        

#endif // G4ANALYSER_HH 
