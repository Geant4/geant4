#include "G4Analyser.hh"
#include <math.h>

G4Analyser::G4Analyser()
  :verboseLevel(2)  {

if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::G4Analyser" << G4endl;
  }

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

void G4Analyser::setInelCsec(G4double csec, 
			     G4bool withn) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::setInelCsec" << G4endl;
  }

  inel_csec = csec; // mb
  withNuclei = withn;

  if (verboseLevel > 1) {
    G4cout << " total inelastic " << inel_csec << G4endl;
  }
}

 
void G4Analyser::setWatchers(const vector<G4NuclWatcher>& watchers) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::setWatchers" << G4endl;
  }

  ana_watchers = watchers;

  G4cout << " watchers set " << watchers.size() << G4endl;
};

void G4Analyser::try_watchers(G4double a, 
			      G4double z, 
			      G4bool if_nucl) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::try_watchers" << G4endl;
  }

  for(G4int iw = 0; iw < ana_watchers.size(); iw++) { 

    if(if_nucl) {
      if(ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a, z); 
    }
    else {
      if(!ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a, z); 
    }; 
  };
};


void G4Analyser::analyse(const G4CollisionOutput& output) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::analyse" << G4endl;
  }

  if(withNuclei) {
    vector<G4InuclNuclei> nucleus = output.getNucleiFragments();

    if(nucleus.size() >= 0) {
      G4int nbig = 0;
      averageNucleiFragments += nucleus.size();

      for(G4int in = 0; in < nucleus.size(); in++) {
	averageExitationEnergy += nucleus[in].getExitationEnergy();

	G4double a = nucleus[in].getA();
	G4double z = nucleus[in].getZ();

	if(in == 0) { 
	  averageA += a; 
	  averageZ += z; 
	};
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
	    zp = 1.0;
	    ap = 1.0;
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

void G4Analyser::printResultsSimple() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::printResultsSimple" << G4endl;
  }

  G4cout << " Number of events " << int(eventNumber + 0.1) << G4endl
	 << " average multiplicity " << averageMultiplicity / eventNumber << G4endl
	 << " average proton number " << averageProtonNumber / eventNumber << G4endl
	 << " average neutron number " << averageNeutronNumber / eventNumber << G4endl
	 << " average nucleon Ekin " << averageNucleonKinEnergy /
    (averageProtonNumber + averageNeutronNumber) << G4endl
	 << " average proton Ekin " << averageProtonKinEnergy / (averageProtonNumber +
								 1.0e-10) << G4endl
	 << " average neutron Ekin " << averageNeutronKinEnergy / (averageNeutronNumber +
								   1.0e-10) << G4endl
	 << " average pion number " << averagePionNumber / eventNumber << G4endl
	 << " average pion Ekin " << averagePionKinEnergy / (averagePionNumber +
							     1.0e-10) << G4endl;
  if(withNuclei) {
    G4cout		   
      << " average Exitation Energy " << 
      averageExitationEnergy / averageNucleiFragments << G4endl
      << " average num of fragments " << averageNucleiFragments / eventNumber << G4endl;
    G4cout << " fission prob. " << fissy_prob / eventNumber << " c.sec " <<
      inel_csec * fissy_prob / eventNumber << G4endl;
  };
};

void G4Analyser::printResults() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::printResults" << G4endl;
  }

  G4cout << " Number of events " << int(eventNumber + 0.1) << G4endl
	 << " average multiplicity " << averageMultiplicity / eventNumber << G4endl
	 << " average proton number " << averageProtonNumber / eventNumber << G4endl
	 << " average neutron number " << averageNeutronNumber / eventNumber << G4endl
	 << " average nucleon Ekin " << averageNucleonKinEnergy /
    (averageProtonNumber + averageNeutronNumber) << G4endl
	 << " average proton Ekin " << averageProtonKinEnergy / (averageProtonNumber +
								 1.0e-10) << G4endl
	 << " average neutron Ekin " << averageNeutronKinEnergy / (averageNeutronNumber +
								   1.0e-10) << G4endl
	 << " average pion number " << averagePionNumber / eventNumber << G4endl
	 << " average pion Ekin " << averagePionKinEnergy / (averagePionNumber +
							     1.0e-10) << G4endl
	 << " average pi+ " << averagePionPl / eventNumber << G4endl
	 << " average pi- " << averagePionMin / eventNumber << G4endl
	 << " average pi0 " << averagePion0 / eventNumber << G4endl;
     		   
  if(withNuclei) {

    G4cout
      << " average A " << averageA / eventNumber << G4endl 		   
      << " average Z " << averageZ / eventNumber << G4endl 		   
      << " average Exitation Energy " << 
      averageExitationEnergy / averageNucleiFragments << G4endl
      << " average num of fragments " << averageNucleiFragments / eventNumber << G4endl;
    G4cout << " fission prob. " << fissy_prob / eventNumber << " c.sec " <<
      inel_csec * fissy_prob / eventNumber << G4endl;
    handleWatcherStatistics();
  };
};

void G4Analyser::handleWatcherStatistics() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::handleWatcherStatistics" << G4endl;
  }

  const G4double small = 1.0e-10;

  if (verboseLevel > 1) {
    G4cout << " ===================================================== " << G4endl;
    G4cout << " ********** Izotop analysis ****************** " << G4endl;
  };

  G4double fgr = 0.0;
  G4double averat = 0.0;
  G4double ave_err = 0.0;
  G4double gl_chsq = 0.0;
  G4double tot_exper = 0.0;
  G4double tot_exper_err = 0.0;
  G4double tot_inucl = 0.0;
  G4double tot_inucl_err = 0.0;
  G4double checked = 0.0;

  for(G4int iw = 0; iw < ana_watchers.size(); iw++) {
    ana_watchers[iw].setInuclCs(inel_csec, eventNumber);
    ana_watchers[iw].print();

    if(ana_watchers[iw].to_check()) {
      pair<G4double, G4double> rat_err = ana_watchers[iw].getAverageRatio();

      averat += rat_err.first;
      ave_err += rat_err.second;
      gl_chsq += ana_watchers[iw].getChsq();   

      pair<G4double, G4double> cs_err = ana_watchers[iw].getExpCs();

      tot_exper += cs_err.first;
      tot_exper_err += cs_err.second;

      pair<G4double, G4double> inucl_cs_err = ana_watchers[iw].getInuclCs();

      tot_inucl += inucl_cs_err.first;
      tot_inucl_err += inucl_cs_err.second;

      G4double iz_checked = ana_watchers[iw].getNmatched();

      if(iz_checked > 0.0) {
	fgr += ana_watchers[iw].getLhood();
	checked += iz_checked;    
      };
    };
  };

  if(checked > 0.0) {
    gl_chsq = sqrt(gl_chsq) / checked;
    averat /= checked;
    ave_err /= checked;
    fgr = pow(10.0, sqrt(fgr / checked)); 
  };

  if (verboseLevel > 1) {
    G4cout << " ===================================================== " << G4endl;
    G4cout << " total exper c.s. " << tot_exper << " err " << tot_exper_err <<
      " tot inucl c.s. " << tot_inucl << " err " << tot_inucl_err << G4endl;
    G4cout << " checked total " << checked << " lhood " << fgr << G4endl
	   << " average ratio " << averat << " err " << ave_err << G4endl
	   << " global chsq " << gl_chsq << G4endl;
  }
}

void G4Analyser::printResultsNtuple() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::printResultsNtuple" << G4endl;
  }

  // Create one line of ACII data. 
  // Several runs should create ntuple for data-analysis 
  G4cout << int(eventNumber + 0.1) << " " <<
    averageMultiplicity / eventNumber << " " << 
    averageProtonNumber / eventNumber << " " <<
    averageNeutronNumber / eventNumber << " " <<
    averageNucleonKinEnergy / (averageProtonNumber + averageNeutronNumber) << " " <<
    averageProtonKinEnergy / (averageProtonNumber + 1.0e-10) << " " <<
    averageNeutronKinEnergy / (averageNeutronNumber + 1.0e-10) << " " <<
    averagePionNumber / eventNumber << " " <<
    averagePionKinEnergy / (averagePionNumber + 1.0e-10) << G4endl;
};


