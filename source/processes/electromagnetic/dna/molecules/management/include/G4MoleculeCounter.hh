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
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4MoleculeCounter_h
#define G4MoleculeCounter_h

#include "G4Molecule.hh"
#include <map>
#include <memory>
#include <set>

struct compDoubleWithPrecision
{
  bool operator()(const double& a, const double& b) const
  {
    if (std::fabs(a - b) < fPrecision)
    {
      return false;
    }
    else
    {
      return a < b;
    }
  }

  static G4ThreadLocal double fPrecision;
};

typedef std::map<G4double, G4int, compDoubleWithPrecision> NbMoleculeAgainstTime;

#if __cplusplus > 199711L && !defined __clang__
#define stdunique_ptr std::unique_ptr
#else
#define stdunique_ptr std::auto_ptr
#endif

typedef stdunique_ptr<std::set<G4double> > RecordedTimes;
typedef std::set<G4double>::iterator RecordedTimesIterator;

class G4MoleculeCounter
{
public:
  typedef std::map<G4Molecule, NbMoleculeAgainstTime> CounterMapType;
  typedef stdunique_ptr<std::vector<G4Molecule> > RecordedMolecules;

  /*
   #if __cplusplus > 199711L && !defined __clang__
   typedef std::unique_ptr<std::vector<G4Molecule> > RecordedMolecules;
   #else
   typedef std::auto_ptr<std::vector<G4Molecule> > RecordedMolecules;
   #endif
   */

protected:
  G4MoleculeCounter();
  virtual ~G4MoleculeCounter()
  { ;}
  static G4ThreadLocal G4MoleculeCounter* fpInstance;

  CounterMapType fCounterMap;
  std::map<const G4MoleculeDefinition*, G4bool> fDontRegister;
  static G4bool fUse;

  G4int fVerbose;

  struct Search
  {
    Search()
    {
      fLowerBoundSet = false;
    }
    CounterMapType::iterator fLastMoleculeSearched;
    NbMoleculeAgainstTime::iterator fLowerBoundTime;
    bool fLowerBoundSet;
  };

  stdunique_ptr<Search> fpLastSearch;

#ifdef MOLECULE_COUNTER_TESTING
public:
#else
protected:
#endif

  friend class G4Molecule;
  virtual void AddAMoleculeAtTime(const G4Molecule&, G4double);
  virtual void RemoveAMoleculeAtTime(const G4Molecule&, G4double);

public:
  static void DeleteInstance();
  static G4MoleculeCounter* Instance();
  static G4MoleculeCounter* GetMoleculeCounter();
  void Initialize();
  static void InitializeInstance();

  G4bool SearchTimeMap(const G4Molecule &molecule);
  int SearchUpperBoundTime(double time, bool sameTypeOfMolecule);

  int GetNMoleculesAtTime(const G4Molecule &molecule, double time);
  inline const NbMoleculeAgainstTime&
  GetNbMoleculeAgainstTime(const G4Molecule &molecule);

  RecordedMolecules GetRecordedMolecules();
  RecordedTimes GetRecordedTimes();

  /*
   * The dynamics of the given molecule won't be saved into memory.
   */
  inline virtual void DontRegister(const G4MoleculeDefinition*);
  inline virtual void RegisterAll();

  /*
   * If the molecule counter is used, it will be called
   * at every creation/deletion of a molecule to
   * to increase/decrease the number at a given time.
   */
  void Use(G4bool flag = true)
  {
    fUse=flag;
  }
  G4bool InUse()
  {
    return fUse;
  }

  inline void SetVerbose(G4int);
  inline G4int GetVerbose();

  /*
   * It sets the min time difference in between two time slices.
   */
  void SetTimeSlice(double);

  virtual void ResetCounter();
};

inline void G4MoleculeCounter::ResetCounter()
{
  fCounterMap.clear();
}

inline const NbMoleculeAgainstTime&
G4MoleculeCounter::GetNbMoleculeAgainstTime(const G4Molecule& molecule)
{
  return fCounterMap[molecule];
}

inline void G4MoleculeCounter::SetVerbose(G4int level)
{
  fVerbose = level;
}

inline G4int G4MoleculeCounter::GetVerbose()
{
  return fVerbose;
}

inline void G4MoleculeCounter::DontRegister(const G4MoleculeDefinition* molDef)
{
  fDontRegister[molDef] = true;
}

inline void G4MoleculeCounter::RegisterAll()
{
  fDontRegister.clear();
}

#endif
