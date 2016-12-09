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
// Author: Mathieu Karamitros

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

#include "G4VMoleculeCounter.hh"
#include <map>
#include <memory>
#include <set>
#include <vector>

//------------------------------------------------------------------------------

struct compDoubleWithPrecision{
  bool operator()(const double& a, const double& b) const{
    if (std::fabs(a - b) < fPrecision){
      return false;
    }
    else{
      return a < b;
    }
  }

  static G4ThreadLocal double fPrecision;
};

//------------------------------------------------------------------------------

typedef std::map<G4double, G4int, compDoubleWithPrecision>
  NbMoleculeAgainstTime;

typedef std::unique_ptr<std::set<G4double> > RecordedTimes;
typedef std::set<G4double>::iterator RecordedTimesIterator;

//------------------------------------------------------------------------------

class G4MoleculeCounter : public G4VMoleculeCounter
{
  //----------------------------------------------------------------------------
public:
  typedef std::map<G4MolecularConfiguration*,
                   NbMoleculeAgainstTime> CounterMapType;
  typedef std::unique_ptr<std::vector<G4MolecularConfiguration*> >
    RecordedMolecules;

  //----------------------------------------------------------------------------
protected:
  G4MoleculeCounter();
  virtual ~G4MoleculeCounter();

  CounterMapType fCounterMap;
  std::map<const G4MoleculeDefinition*, G4bool> fDontRegister;

  G4int fVerbose;
  G4bool fCheckTimeIsConsistentWithScheduler;

  struct Search
  {
    Search(){
      fLowerBoundSet = false;
    }
    CounterMapType::iterator fLastMoleculeSearched;
    NbMoleculeAgainstTime::iterator fLowerBoundTime;
    bool fLowerBoundSet;
  };

  std::unique_ptr<Search> fpLastSearch;

#ifdef MOLECULE_COUNTER_TESTING
public:
#else
protected:
#endif

  friend class G4Molecule;
  friend class G4VMoleculeCounter;

  //----------------------------------------------------------------------------
public:
  static G4MoleculeCounter* Instance();
  
  //----------------------------------------------------------------------------
public:
  void Initialize() override;
  
  void ResetCounter() override;
  
  void AddAMoleculeAtTime(G4MolecularConfiguration*,
                          G4double time,
                          const G4ThreeVector* position = nullptr,
                          int number = 1) override;
  void RemoveAMoleculeAtTime(G4MolecularConfiguration*,
                             G4double time,
                             const G4ThreeVector* position = nullptr,
                             int number = 1) override;
  
  /* The dynamics of the given molecule won't be saved into memory.*/
  inline void DontRegister(const G4MoleculeDefinition*) override;
  inline bool IsRegistered(const G4MoleculeDefinition*) override;
  inline void RegisterAll() override;

  //----------------------------------------------------------------------------
public:
  int GetNMoleculesAtTime(G4MolecularConfiguration* molecule,
                          double time);
  inline const NbMoleculeAgainstTime&
  GetNbMoleculeAgainstTime(G4MolecularConfiguration* molecule);

  RecordedMolecules GetRecordedMolecules();
  RecordedTimes GetRecordedTimes();

  inline void SetVerbose(G4int);
  inline G4int GetVerbose();

  /*
   * It sets the min time difference in between two time slices.
   */
  void SetTimeSlice(double);

  
  void Dump();

  inline G4bool IsTimeCheckedForConsistency() const
  {
    return fCheckTimeIsConsistentWithScheduler;
  }

  inline void CheckTimeForConsistency(G4bool flag)
  {
    fCheckTimeIsConsistentWithScheduler = flag;
  }
  
  //----------------------------------------------------------------------------
protected:
  G4bool SearchTimeMap(G4MolecularConfiguration* molecule);
  int SearchUpperBoundTime(double time, bool sameTypeOfMolecule);
};

//------------------------------------------------------------------------------

inline void G4MoleculeCounter::ResetCounter()
{
  if(fVerbose)
  {
      G4cout << " ---> G4MoleculeCounter::ResetCounter" << G4endl;
  }
  fCounterMap.clear();
  fpLastSearch.reset(0);
}

//------------------------------------------------------------------------------

inline const NbMoleculeAgainstTime&
G4MoleculeCounter::GetNbMoleculeAgainstTime(G4MolecularConfiguration* molecule)
{
  return fCounterMap[molecule];
}

//------------------------------------------------------------------------------

inline void G4MoleculeCounter::SetVerbose(G4int level)
{
  fVerbose = level;
}

//------------------------------------------------------------------------------

inline G4int G4MoleculeCounter::GetVerbose()
{
  return fVerbose;
}

//------------------------------------------------------------------------------

inline void
G4MoleculeCounter::DontRegister(const G4MoleculeDefinition* molDef)
{
  fDontRegister[molDef] = true;
}

//------------------------------------------------------------------------------

inline bool
G4MoleculeCounter::IsRegistered(const G4MoleculeDefinition* molDef)
{
  if(fDontRegister.find(molDef) == fDontRegister.end()) return true;
  return false;
}

//------------------------------------------------------------------------------

inline void G4MoleculeCounter::RegisterAll()
{
  fDontRegister.clear();
}

#endif
