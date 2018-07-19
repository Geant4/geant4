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
/// \file parallel/ThreadsafeScorers/include/TSRun.hh
/// \brief Definition of the TSRun class
//
//
// $Id: TSRun.hh 93110 2015-11-05 08:37:42Z jmadsen $
//
//
/// TSRun contains three collections of hits maps: a thread-local hits map,
///     a global atomic hits map (implemented as a static since TSRun is
///     implemented as a thread-local instance), and a global "mutex" hits map
///     (also implemented as a static). The thread-local hits map is the
///     same as you will find in many other examples. The atomics hits map
///     is the purpose of this example. Code-wise, the implementation looks
///     extremely similar to the thread-local version with the 3 primary
///     exceptions: (1) construction - there should only be one instance so
///     it should be a static member variable or a pointer/reference to a
///     single instance elsewhere in the code (stored in ActionInitialization,
///     for instance); (2) It does not need to, nor should be, summed in
///     G4Run::Merge(); and (3) destruction -- it should only be cleared by
///     the master thread since there is only one instance.
/// A "mutex" hits map is also included as reference for checking the results
///     accumulated by the thread-local hits maps and atomic hits maps. The
///     differences w.r.t. this hits maps are computed in
///     TSRunAction::EndOfRunAction
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef tsrun_h
#define tsrun_h 1

#include "globals.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4TAtomicHitsMap.hh"

#include <vector>

class G4Event;

class TSRun : public G4Run
{
public:
    typedef std::map<G4int, G4double> MutexHitsMap_t;

public:
    TSRun(const G4String&);
    virtual ~TSRun();

    // virtual method from G4Run.
    // The method is overriden in this class for scoring.
    virtual void RecordEvent(const G4Event*);

    // Access methods for scoring information.
    // - Get HitsMap of this RUN.
    G4THitsMap<G4double>* GetHitsMap(const G4String& collname) const;
    G4TAtomicHitsMap<G4double>* GetAtomicHitsMap(const G4String&) const;
    MutexHitsMap_t* GetMutexHitsMap(const G4String&) const;

    void ConstructMFD(const G4String&);

    virtual void Merge(const G4Run*);

private:
    std::vector<G4String> fCollNames;
    std::vector<G4int> fCollIDs;
    std::vector<G4THitsMap<G4double>*> fRunMaps;
    static std::vector<G4TAtomicHitsMap<G4double>*> fAtomicRunMaps;
    static std::map<G4String, MutexHitsMap_t> fMutexRunMaps;

};

#endif
