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
// $Id: G4CascadeHistory.hh 67796 2013-03-08 06:18:39Z mkelsey $
//
// G4CascadeHistory: Container to record all particles produced during
// cascade, with daughters; printout is formatted hierarchically.

#ifndef G4CASCADE_HISTORY_HH
#define G4CASCADE_HISTORY_HH

#include "globals.hh"
#include "G4CascadParticle.hh"
#include <iosfwd>
#include <set>
#include <vector>


class G4CascadeHistory {
public:
  G4CascadeHistory(G4int verbose=0) : verboseLevel(verbose) {;}
  ~G4CascadeHistory() {;}		// *** Do not want subclasses ***

  void setVerboseLevel(G4int verbose=0) { verboseLevel = verbose; }

  // Reset buffers for new event
  void Clear();

  // Add particle to history list, assigning ID number (non-const input)
  G4int AddEntry(G4CascadParticle& cpart);

  // Record full interaction vertex (non-const input)
  G4int AddVertex(G4CascadParticle& cpart, std::vector<G4CascadParticle>& daug);

  // Discard particle reabsorbed during cascade
  void DropEntry(const G4CascadParticle& cpart);

  // Report cascade structure hierarchically
  void Print(std::ostream& os) const;

protected:
  struct HistoryEntry {
    G4CascadParticle cpart;
    G4int n;
    G4int dId[10];		// Must be fixed size for allocation in vector

    HistoryEntry() { clear(); };
    HistoryEntry(const G4CascadParticle& cp) : cpart(cp) { clear(); }
    void clear();
  };

  // Assign ID number to particle (non-const input)
  void AssignHistoryID(G4CascadParticle& cpart);

  // Populate list of daughters in interaction, adding them to history list
  void FillDaughters(G4int iEntry, std::vector<G4CascadParticle>& daug);

  // Add single-line report for particle, along with daughters
  void PrintEntry(std::ostream& os, G4int iEntry) const;
  void PrintParticle(std::ostream& os, const G4CascadParticle& cpart) const;
  G4bool PrintingDone(G4int iEntry) const {
    return (entryPrinted.find(iEntry) != entryPrinted.end());
  }

  // Derive target of cascade step from particle and daughters
  const char* GuessTarget(const HistoryEntry& entry) const;

  G4int size() const { return (G4int)theHistory.size(); }

private:
  G4int verboseLevel;

  std::vector<HistoryEntry> theHistory;		// List of particles and daughters
  mutable std::set<G4int> entryPrinted;		// Particle indices already reported

  // No copying allowed
  G4CascadeHistory(const G4CascadeHistory& rhs);
  G4CascadeHistory& operator=(const G4CascadeHistory& rhs);
};

// Inline reporting
std::ostream& operator<<(std::ostream& os, const G4CascadeHistory& history);

#endif	/* G4CASCADE_HISTORY_HH */
