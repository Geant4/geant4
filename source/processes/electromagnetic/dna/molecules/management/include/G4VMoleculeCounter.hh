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
//  G4VMoleculeCounter.hh
//  Geant4
//
//  Created by Mathieu Karamitros on 02/11/2016.
//  Modified by Christian Velten on 10/27/2024.
//
//
#ifndef G4VMOLECULECOUNTER_HH
#define G4VMOLECULECOUNTER_HH 1

#include "G4MoleculeCounterTimeComparer.hh"
#include "G4VMoleculeCounterInternalBase.hh"

#include <map>
#include <memory>

class G4MolecularConfiguration;
class G4MoleculeDefinition;
class G4Track;
class G4StepPoint;

//------------------------------------------------------------------------------

class G4VMoleculeCounter : public G4VMoleculeCounterInternalBase
{
    template<typename>
    friend class G4VUserMoleculeCounter;

  public:
    enum MoleculeCounterType
    {
      Other,
      Basic,
      Mesoscopic,
    };
    struct G4VMoleculeCounterIndex : public G4VMoleculeCounterIndexInterface
    {
        virtual ~G4VMoleculeCounterIndex() = default;
        virtual G4bool operator<(G4VMoleculeCounterIndex const&) const = 0;
        virtual G4bool operator==(G4VMoleculeCounterIndex const&) const = 0;
        virtual G4String GetInfo() const = 0;
        virtual const G4MolecularConfiguration* GetMolecule() const = 0;
    };

  private:
    G4VMoleculeCounter();
    G4VMoleculeCounter(const G4String&, MoleculeCounterType = MoleculeCounterType::Other);
    G4VMoleculeCounter(G4VMoleculeCounter const&) = delete;
    void operator=(G4VMoleculeCounter const& x) = delete;

  public:
    ~G4VMoleculeCounter() override = default;

  public:
    virtual std::unique_ptr<G4VMoleculeCounterIndex> BuildIndex(const G4Track*) const = 0;
    virtual std::unique_ptr<G4VMoleculeCounterIndex> BuildIndex(const G4Track*,
                                                                const G4StepPoint*) const = 0;
    virtual std::unique_ptr<G4VMoleculeCounterIndex> BuildSimpleIndex(const G4MolecularConfiguration*) const = 0;

    virtual void AddMolecule(std::unique_ptr<G4VMoleculeCounterIndex>, G4double, G4int = 1) = 0;
    virtual void RemoveMolecule(std::unique_ptr<G4VMoleculeCounterIndex>, G4double, G4int = 1) = 0;

    virtual std::set<const G4MolecularConfiguration*> GetRecordedMolecules() const = 0;
    std::set<G4double> GetRecordedTimes() const override = 0;

    virtual void SchedulerFinalizedTracking() = 0;

  protected:
    MoleculeCounterType fType{MoleculeCounterType::Other};

    G4bool fSensitiveToStepping{false};
    G4bool fNegativeCountsAreFatal{true};

    std::set<const G4MoleculeDefinition*> fIgnoredMolecules{};
    std::set<const G4MolecularConfiguration*> fIgnoredReactants{};

  public:
    MoleculeCounterType GetType() const;

    G4bool GetSensitiveToStepping() const;
    void SetSensitiveToStepping(G4bool = true);

    G4bool GetNegativeCountsAreFatal() const;

    // Ignore certain molecules from counting
    void IgnoreMolecule(const G4MoleculeDefinition*);
    void IgnoreReactant(const G4MolecularConfiguration*);
    void RegisterAll();

    G4bool IsReactantIgnored(const G4MoleculeDefinition*) const;
    G4bool IsReactantIgnored(const G4MolecularConfiguration*) const;

  protected:
    void SetNegativeCountsAreFatal(G4bool);
};

//------------------------------------------------------------------------------

inline G4VMoleculeCounter::MoleculeCounterType G4VMoleculeCounter::GetType() const
{
  return fType;
}

inline G4bool G4VMoleculeCounter::GetSensitiveToStepping() const
{
  return fSensitiveToStepping;
}

inline G4bool G4VMoleculeCounter::GetNegativeCountsAreFatal() const
{
  return fNegativeCountsAreFatal;
}
inline void G4VMoleculeCounter::SetNegativeCountsAreFatal(G4bool flag)
{
  fNegativeCountsAreFatal = flag;
}

inline void G4VMoleculeCounter::IgnoreMolecule(const G4MoleculeDefinition* molecule)
{
  fIgnoredMolecules.insert(molecule);
}

inline void G4VMoleculeCounter::IgnoreReactant(const G4MolecularConfiguration* reactant)
{
  fIgnoredReactants.insert(reactant);
}

inline void G4VMoleculeCounter::RegisterAll()
{
  fIgnoredMolecules.clear();
  fIgnoredReactants.clear();
}

#endif
