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
// $Id: G4CascadeHistory.cc 67796 2013-03-08 06:18:39Z mkelsey $
//
// G4CascadeHistory: Container to record all particles produced during
// cascade, with daughters; printout is formatted hierarchically.

#include "G4CascadeHistory.hh"
#include "G4CascadParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclParticle.hh"
#include <iostream>
#include <iomanip>
#include <ios>


// Constructors for individual entries (vertices)

void G4CascadeHistory::HistoryEntry::clear() {
  for (G4int i=0; i<10; i++) dId[i] = -1;
  n = 0;
}


// Initialize buffers for new event

void G4CascadeHistory::Clear() {
  if (verboseLevel>1) G4cout << " >>> G4CascadeHistory::Clear" << G4endl;
  theHistory.clear();
  entryPrinted.clear();
}


// Record a new cascade vertex (particle and daughters)

G4int G4CascadeHistory::AddVertex(G4CascadParticle& cpart,
				  std::vector<G4CascadParticle>& daug) {
  if (verboseLevel>1) G4cout << " >>> G4CascadeHistory::AddVertex" << G4endl;

  // Create new entry for vertex or update particle kinematics
  G4int id = AddEntry(cpart);
  FillDaughters(id, daug);

  if (verboseLevel>3) {
    G4cout << " entry " << id << " " << &theHistory[id] << " got "
	   << theHistory[id].n << " daughters:";
    for (G4int i=0; i<theHistory[id].n; i++) {
      G4cout << " " << theHistory[id].dId[i];
    }
    G4cout << G4endl;
  }

  return id;
}

void 
G4CascadeHistory::FillDaughters(G4int iEntry, std::vector<G4CascadParticle>& daug) {
  G4int nDaug = (G4int)daug.size();

  if (verboseLevel>1) 
    G4cout << " >>> G4CascadeHistory::FillDaughters " << iEntry << G4endl;

  // NOTE:  Cannot use reference to element, as push_back can invalidate refs!
  theHistory[iEntry].clear();

  theHistory[iEntry].n = nDaug;
  for (G4int i=0; i<nDaug; i++) {
    G4int id = AddEntry(daug[i]);
    theHistory[iEntry].dId[i] = id;
  }

  if (verboseLevel>3) {
    G4cout << " got " << theHistory[iEntry].n << " daughters:";
    for (G4int i=0; i<theHistory[iEntry].n; i++) {
      G4cout << " " << theHistory[iEntry].dId[i];
    }
    G4cout << G4endl;
  }
}

// Add particle to history list, or update if already recorded

G4int G4CascadeHistory::AddEntry(G4CascadParticle& cpart) {
  AssignHistoryID(cpart);		// Make sure particle has index

  G4int id = cpart.getHistoryId();
  if (id < size()) {
    if (verboseLevel>2) 
      G4cout << " AddEntry updating " << id << " " << &theHistory[id] << G4endl;
    theHistory[id].cpart = cpart;		// Copies kinematics
  } else {
    theHistory.push_back(HistoryEntry(cpart));
    if (verboseLevel>2) 
      G4cout << " AddEntry creating " << id << " " << &theHistory.back() << G4endl;
  }

  if (verboseLevel>3) G4cout << theHistory[id].cpart << G4endl;	// Sanity check

  return id;
}

// Discard particle reabsorbed during cascade

void G4CascadeHistory::DropEntry(const G4CascadParticle& cpart) {
  if (verboseLevel>1) G4cout << " >>> G4CascadeHistory::DropEntry" << G4endl;

  
  G4int id = cpart.getHistoryId();	// Particle must appear in history
  if (id>=0) theHistory[id].n = -1;	// Special flag for absorbed particle
}

// Check if particle already in history, assign ID if not there

void G4CascadeHistory::AssignHistoryID(G4CascadParticle& cpart) {
  if (cpart.getHistoryId() >= 0) return;		// ID already assigned

  if (verboseLevel>2) {
    G4cout << " >>> G4CascadeHistory::NewHistoryID assigning ID "
	   << size() << G4endl;
  }

  cpart.setHistoryId(size());
}


// Generate hierarchical (indented) report of history

std::ostream& operator<<(std::ostream& os, const G4CascadeHistory& history) {
  history.Print(os);
  return os;
}

void G4CascadeHistory::Print(std::ostream& os) const {
  if (verboseLevel) os << " >>> G4CascadeHistory::Print" << std::endl;

  os << " Cascade structure: vertices, (-O-) exciton, (***) outgoing"
     << std::endl;

  for (G4int i=0; i<size(); i++) {
    if (!PrintingDone(i)) PrintEntry(os, i);    
  }
}

// Add single-line report for particle, along with daughters
  
void G4CascadeHistory::PrintEntry(std::ostream& os, G4int iEntry) const {
  if (iEntry >= size()) return;			// Skip nonexistent entry
  if (PrintingDone(iEntry)) return;		// Skip entry already reported

  entryPrinted.insert(iEntry);

  const HistoryEntry& entry = theHistory[iEntry];	// For convenience
  const G4CascadParticle& cpart = entry.cpart;

  G4int indent = cpart.getGeneration()*2;

  // Index and indentation of cascade vertex
  std::ios::fmtflags osFlags = os.flags();
  os.setf(std::ios::left);	// Pushes all blanks to right end of output
  os << "#" << std::setw(3+indent) << iEntry;
  os.flags(osFlags);

  os << cpart.getParticle().getDefinition()->GetParticleName()
     << " p " << cpart.getMomentum() << " (cosTh "
     << cpart.getMomentum().vect().unit().z() << ")"
     << " @ " << cpart.getPosition()
     << " zone " << cpart.getCurrentZone();

  // Flag as final-state particle or report daughters iteratively
  os << " (" << GuessTarget(entry) << ")";
  if (entry.n > 0) {
    os << " -> N=" << entry.n << std::endl;
    for (G4int i=0; i<entry.n; i++) {
      PrintEntry(os, entry.dId[i]);
    }
  } else os << std::endl;
}

// Derive target of cascade step from particle and daughters

const char* G4CascadeHistory::
GuessTarget(const G4CascadeHistory::HistoryEntry& entry) const {
  if (verboseLevel>2) G4cout << " >>> G4CascadeHistory::GuessTarget" << G4endl;

  if (entry.n < 0) return "-O-";	// Exciton or trapped-decay
  if (entry.n == 0) return "***";	// Outgoing (final state) particle

  const G4CascadParticle& cpart = entry.cpart;	// For convenience
  if (verboseLevel>3) G4cout << "cpart: " << cpart;

  // Compute baryon number and charge from daughters minus projectile
  G4int targetB = -cpart.getParticle().baryon();
  G4int targetQ = (G4int)-cpart.getParticle().getCharge();

  for (G4int i=0; i<entry.n; i++) {
    const G4CascadParticle& cdaug = theHistory[entry.dId[i]].cpart;
    if (verboseLevel>3) 
      G4cout << "cdaug " << i << " ID " << entry.dId[i] << ": " << cdaug;

    targetB += cdaug.getParticle().baryon();
    targetQ += (G4int)cdaug.getParticle().getCharge();
  }

  // Target possibilities are proton, neutron or dibaryon (pp, nn, pn)
  if (targetB==1 && targetQ==0) return "n";
  if (targetB==1 && targetQ==1) return "p";
  if (targetB==2 && targetQ==0) return "nn";
  if (targetB==2 && targetQ==1) return "pn";
  if (targetB==2 && targetQ==2) return "pp";

  if (verboseLevel>2) {
    G4cout << " ERROR identifying target: deltaB " << targetB
	   << " deltaQ " << targetQ << " from\n" << cpart << " to" << G4endl;
    for (G4int j=0; j<entry.n; j++) {
      G4cout << theHistory[entry.dId[j]].cpart;
    }
  }

  return "BAD TARGET";		// Should not get here if EPCollider worked right
}
