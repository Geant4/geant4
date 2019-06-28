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

#pragma once

#include "G4VMoleculeCounter.hh"
#include <map>
#include <memory>
#include <set>
#include <vector>

//------------------------------------------------------------------------------

namespace G4 {
namespace MoleculeCounter {
struct TimePrecision
{
    bool operator()(const double& a, const double& b) const;
    static G4ThreadLocal double fPrecision;
};
}
}
using NbMoleculeAgainstTime = std::map<G4double, G4int, G4::MoleculeCounter::TimePrecision>;
using RecordedTimes = std::unique_ptr<std::set<G4double>>;

//------------------------------------------------------------------------------

class G4MoleculeCounter : public G4VMoleculeCounter
{
    //----------------------------------------------------------------------------
public:
    using ReactantList = std::vector<Reactant*>;
    using CounterMapType = std::map<Reactant*, NbMoleculeAgainstTime>;
    using RecordedMolecules = std::unique_ptr<ReactantList>;

    static G4MoleculeCounter* Instance();

    void Initialize() override;
    void ResetCounter() override;

    /* The dynamics of the given molecule won't be saved into memory.*/
    void DontRegister(const G4MoleculeDefinition*) override;
    bool IsRegistered(const G4MoleculeDefinition*) override;
    void RegisterAll() override;

    //----------------------------------------------------------------------------

    int GetNMoleculesAtTime(Reactant* molecule, double time);
    const NbMoleculeAgainstTime& GetNbMoleculeAgainstTime(Reactant* molecule);

    RecordedMolecules GetRecordedMolecules();
    RecordedTimes GetRecordedTimes();

    void SetVerbose(G4int);
    G4int GetVerbose();

    /* It sets the min time difference in between two time slices. */
    static void SetTimeSlice(double);

    void Dump();

    G4bool IsTimeCheckedForConsistency() const;
    void CheckTimeForConsistency(G4bool flag);

#ifdef MOLECULE_COUNTER_TESTING
public:
#else
protected:
#endif
    void AddAMoleculeAtTime(Reactant*,
                            G4double time,
                            const G4ThreeVector* position = nullptr,
                            int number = 1) override;
    void RemoveAMoleculeAtTime(Reactant*,
                               G4double time,
                               const G4ThreeVector* position = nullptr,
                               int number = 1) override;

    //----------------------------------------------------------------------------
protected:
    G4bool SearchTimeMap(Reactant* molecule);
    int SearchUpperBoundTime(double time, bool sameTypeOfMolecule);

protected:
    G4MoleculeCounter();
    ~G4MoleculeCounter() override;

    CounterMapType fCounterMap;
    std::map<const G4MoleculeDefinition*, G4bool> fDontRegister;

    G4int fVerbose;
    G4bool fCheckTimeIsConsistentWithScheduler;

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

    std::unique_ptr<Search> fpLastSearch;

    friend class G4Molecule;
    friend class G4VMoleculeCounter;
};