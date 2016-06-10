//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4Analyser.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100726  M. Kelsey -- Use references for fetched lists
// 20101010  M. Kelsey -- Migrate to integer A and Z
// 20101019  M. Kelsey -- CoVerity report, unitialized constructor

#include "G4Analyser.hh"
#include <cmath>
#include <iomanip>

G4Analyser::G4Analyser()
  : verboseLevel(0), eventNumber(0.0), averageMultiplicity(0.0),
    averageProtonNumber(0.0), averageNeutronNumber(0.0),
    averagePionNumber(0.0),  averageNucleonKinEnergy(0.0),
    averageProtonKinEnergy(0.0), averageNeutronKinEnergy(0.0),
    averagePionKinEnergy(0.0), averageExitationEnergy(0.0),
    averageOutgoingNuclei(0.0), fissy_prob(0.0), averagePionPl(0.0),
    averagePionMin(0.0), averagePion0(0.0), averageA(0.0), averageZ(0.0),
    inel_csec(0.0), withNuclei(false) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::G4Analyser" << G4endl;
  }
}

void G4Analyser::setInelCsec(G4double csec, G4bool withn) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::setInelCsec" << G4endl;
  }

  inel_csec = csec; // mb
  withNuclei = withn;

  if (verboseLevel > 3) {
    G4cout << " total inelastic " << inel_csec << G4endl;
  }
}

void G4Analyser::setWatchers(const std::vector<G4NuclWatcher>& watchers) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::setWatchers" << G4endl;
  }

  ana_watchers = watchers;
  if (verboseLevel > 3) {
    G4cout << " watchers set " << watchers.size() << G4endl;
  }
}

void G4Analyser::try_watchers(G4int a, G4int z, G4bool if_nucl) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::try_watchers" << G4endl;
  }

  for (G4int iw = 0; iw < G4int(ana_watchers.size()); iw++) { 

    if (if_nucl) {

      if (ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a, z); 

    } else {

      if (!ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a, z); 
    }
  }
}


void G4Analyser::analyse(const G4CollisionOutput& output) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::analyse" << G4endl;
  }

  if (withNuclei) {
    const std::vector<G4InuclNuclei>& nucleus = output.getOutgoingNuclei();

    //    if (nucleus.size() >= 0) {
    if (nucleus.size() > 0) {
      G4int nbig = 0;
      averageOutgoingNuclei += nucleus.size();

      for (G4int in = 0; in < G4int(nucleus.size()); in++) {
	averageExitationEnergy += nucleus[in].getExitationEnergy();

	G4int a = nucleus[in].getA();
	G4int z = nucleus[in].getZ();

	if (in == 0) { 
	  averageA += a; 
	  averageZ += z; 
	};

	if (a > 10) nbig++;
	try_watchers(a, z, true);
      };

      if (nbig > 1) fissy_prob += 1.0;
      eventNumber += 1.0;
      const std::vector<G4InuclElementaryParticle>& particles =
	output.getOutgoingParticles();
      averageMultiplicity += particles.size();

      for (G4int i = 0; i < G4int(particles.size()); i++) {
	G4int ap = 0;
	G4int zp = 0;

	if (particles[i].nucleon()) {
	  averageNucleonKinEnergy += particles[i].getKineticEnergy();

	  if (particles[i].type() == 1) {
	    zp = 1;
	    ap = 1;
	    averageProtonNumber += 1.0;
	    averageProtonKinEnergy += particles[i].getKineticEnergy();

	  } else {
	    ap = 1;
	    zp = 0;
	    averageNeutronNumber += 1.0;
	    averageNeutronKinEnergy += particles[i].getKineticEnergy();
	  };  

	} else if (particles[i].pion()) {
	  averagePionKinEnergy += particles[i].getKineticEnergy();
	  averagePionNumber += 1.0;
	  ap = 0;

	  if (particles[i].type() == 3) {
	    zp = 1;
	    averagePionPl += 1.0;

	  } else if (particles[i].type() == 5) {  
	    zp = -1;
	    averagePionMin += 1.0;

	  } else if (particles[i].type() == 7) { 
	    zp = 0;
	    averagePion0 += 1.0;
	  };
	};
	try_watchers(ap, zp, false);
      };
    };

  } else {
    eventNumber += 1.0;
    const std::vector<G4InuclElementaryParticle>& particles =
      output.getOutgoingParticles();
    averageMultiplicity += particles.size();

    for (G4int i = 0; i < G4int(particles.size()); i++) {

      if (particles[i].nucleon()) {
	averageNucleonKinEnergy += particles[i].getKineticEnergy();

	if (particles[i].type() == 1) {
	  averageProtonNumber += 1.0;
	  averageProtonKinEnergy += particles[i].getKineticEnergy();

	} else {
	  averageNeutronNumber += 1.0;
	  averageNeutronKinEnergy += particles[i].getKineticEnergy();
	}

      } else if (particles[i].pion()) {
	averagePionKinEnergy += particles[i].getKineticEnergy();
	averagePionNumber += 1.0;
      }
    }
  }
}

void G4Analyser::printResultsSimple() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::printResultsSimple" << G4endl;
  }

  G4cout << " Number of events " << G4int(eventNumber + 0.1) << G4endl
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
  if (withNuclei) {
    G4cout		   
      << " average Exitation Energy " << 
      averageExitationEnergy / averageOutgoingNuclei << G4endl
      << " average num of fragments " << averageOutgoingNuclei / eventNumber << G4endl;
    G4cout << " fission prob. " << fissy_prob / eventNumber << " c.sec " <<
      inel_csec * fissy_prob / eventNumber << G4endl;
  }
}

void G4Analyser::printResults() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::printResults" << G4endl;
  }

  G4cout << " Number of events " << G4int(eventNumber + 0.1) << G4endl
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
     		   
  if (withNuclei) {
    G4cout
      << " average A " << averageA / eventNumber << G4endl 		   
      << " average Z " << averageZ / eventNumber << G4endl 		   
      << " average Exitation Energy " << 
      averageExitationEnergy / averageOutgoingNuclei << G4endl
      << " average num of fragments " << averageOutgoingNuclei / eventNumber << G4endl;
    G4cout << " fission prob. " << fissy_prob / eventNumber << " c.sec " <<
      inel_csec * fissy_prob / eventNumber << G4endl;
    handleWatcherStatistics();
  }
}

void G4Analyser::handleWatcherStatistics() {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Analyser::handleWatcherStatistics" << G4endl;
  }

  // const G4double small = 1.0e-10;

  if (verboseLevel > 3) {
    G4cout << " >>>Izotop analysis:" << G4endl;
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

  for (G4int iw = 0; iw < G4int(ana_watchers.size()); iw++) {
    ana_watchers[iw].setInuclCs(inel_csec, G4int(eventNumber));
    ana_watchers[iw].print();

    if (ana_watchers[iw].to_check()) {
      std::pair<G4double, G4double> rat_err = ana_watchers[iw].getAverageRatio();
      averat += rat_err.first;
      ave_err += rat_err.second;
      gl_chsq += ana_watchers[iw].getChsq();   
      std::pair<G4double, G4double> cs_err = ana_watchers[iw].getExpCs();
      tot_exper += cs_err.first;
      tot_exper_err += cs_err.second;
      std::pair<G4double, G4double> inucl_cs_err = ana_watchers[iw].getInuclCs();
      tot_inucl += inucl_cs_err.first;
      tot_inucl_err += inucl_cs_err.second;
      G4double iz_checked = ana_watchers[iw].getNmatched();

      if (iz_checked > 0.0) {
	fgr += ana_watchers[iw].getLhood();
	checked += iz_checked;    
      };
    };
  };

  if (checked > 0.0) {
    gl_chsq = std::sqrt(gl_chsq) / checked;
    averat /= checked;
    ave_err /= checked;
    fgr = std::pow(10.0, std::sqrt(fgr / checked)); 
  };

  if (verboseLevel > 3) {
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
  G4cout <<
    std::setw(15) << int(eventNumber + 0.1) <<
    std::setw(15) << averageMultiplicity / eventNumber << 
    std::setw(15) << averageProtonNumber / eventNumber <<
    std::setw(15) << averageNeutronNumber / eventNumber << " " <<
    std::setw(15) << averageNucleonKinEnergy / (averageProtonNumber + averageNeutronNumber) << " " <<
    std::setw(15) << averageProtonKinEnergy / (averageProtonNumber + 1.0e-10) << " " <<
    std::setw(15) << averageNeutronKinEnergy / (averageNeutronNumber + 1.0e-10) << " " <<
    std::setw(15) << averagePionNumber / eventNumber << " " <<
    std::setw(15) << averagePionKinEnergy / (averagePionNumber + 1.0e-10) << G4endl;
}
