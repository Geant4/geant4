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
// Author: Christian Velten (2025)

#ifndef G4VMoleculeCounterInternalBaseBASE_HH
#define G4VMoleculeCounterInternalBaseBASE_HH 1

#include "G4Exception.hh"
#include "G4MoleculeCounterTimeComparer.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include <iomanip>
#include <map>
#include <memory>
#include <set>

//------------------------------------------------------------------------------
namespace G4
{
namespace MoleculeCounter
{
struct FixedTimeComparer
{
    G4bool operator()(const G4double& a, const G4double& b) const;
    static G4ThreadLocal G4double fPrecision;
};
}  // namespace MoleculeCounter
}  // namespace G4

//------------------------------------------------------------------------------

using InnerCounterMapType = std::map<G4double, G4int, G4MoleculeCounterTimeComparer>;

//------------------------------------------------------------------------------

class G4VMoleculeCounterInternalBase
{
    friend class G4VMoleculeCounter;
    friend class G4VMoleculeReactionCounter;
    // only allow these classes to derive

  public:
    struct G4VMoleculeCounterIndexInterface
    {
        virtual ~G4VMoleculeCounterIndexInterface() = default;
        virtual G4String GetInfo() const = 0;
    };

  private:
    G4VMoleculeCounterInternalBase();
    G4VMoleculeCounterInternalBase(const G4String&);
    G4VMoleculeCounterInternalBase(G4VMoleculeCounterInternalBase const&) = delete;
    void operator=(G4VMoleculeCounterInternalBase const& x) = delete;

  public:
    virtual ~G4VMoleculeCounterInternalBase() = default;

  public:
    virtual void Initialize() = 0;
    virtual void InitializeUser() = 0;
    virtual void ResetCounter() = 0;
    virtual void Dump() const = 0;
    virtual void DumpCounterMapIndices() const = 0;

    virtual std::set<G4double> GetRecordedTimes() const = 0;

    virtual void AbsorbCounter(const G4VMoleculeCounterInternalBase*) = 0;

  private:
    static G4ThreadLocal G4int _createdCounters;

  protected:
    G4bool fIsInitialized{false};

    G4int fId;
    G4int fManagedId{-1};
    G4String fName{};

    G4int fVerbose{0};

    G4double fActiveLowerBound{0}, fActiveUpperBound{std::numeric_limits<G4double>::max()};
    G4bool fActiveLowerBoundInclusive{true}, fActiveUpperBoundInclusive{true};

    G4bool fCheckTimeIsConsistentWithScheduler{true};
    G4bool fCheckRecordedTimesAreConsistent{true};

    G4MoleculeCounterTimeComparer fTimeComparer{};

  public:
    G4int GetId() const;
    void SetManagedId(G4int);
    G4int GetManagedId() const;
    const G4String& GetName() const;

    G4int GetVerbose() const;
    void SetVerbose(G4int);

    // Set counter active w.r.t. time
    G4double GetActiveLowerBound() const;
    void SetActiveLowerBound(G4double, G4bool = true);
    G4double GetActiveUpperBound() const;
    void SetActiveUpperBound(G4double, G4bool = true);
    G4bool GetActiveLowerBoundInclusive() const;
    G4bool GetActiveUpperBoundInclusive() const;
    G4bool IsTimeBelowLowerBound(G4double) const;
    G4bool IsTimeAboveUpperBound(G4double) const;
    G4bool IsActiveAtGlobalTime(G4double) const;

    G4bool GetCheckTimeConsistencyWithScheduler() const; // w.r.t. scheduler time
    void SetCheckTimeConsistencyWithScheduler(G4bool = true);
    G4bool GetCheckRecordedTimeConsistency() const; // w.r.t last recorded time
    void SetCheckRecordedTimeConsistency(G4bool = true);

    const G4MoleculeCounterTimeComparer& GetTimeComparer() const;
    void SetTimeComparer(const G4MoleculeCounterTimeComparer&);

  public:
    static void SetFixedTimePrecision(G4double);
};

//------------------------------------------------------------------------------

inline G4int G4VMoleculeCounterInternalBase::GetId() const
{
  return fId;
}

inline void G4VMoleculeCounterInternalBase::SetManagedId(G4int id)
{
  if (fManagedId > -1) {
    G4ExceptionDescription description;
    description << "Someone is trying to change the managed id of this counter but it was already "
                   "changed from -1!\n";
    description << "  Id: " << fManagedId << "\n";
    description << "Name: " << fName << "\n";
    G4Exception("G4VMoleculeCounterInternalBase::SetManagedId", "MOLCTR000", FatalException, description);
  }
  fManagedId = id;
}
inline G4int G4VMoleculeCounterInternalBase::GetManagedId() const
{
  return fManagedId;
}

inline const G4String& G4VMoleculeCounterInternalBase::GetName() const
{
  return fName;
}

inline G4int G4VMoleculeCounterInternalBase::GetVerbose() const
{
  return fVerbose;
}
inline void G4VMoleculeCounterInternalBase::SetVerbose(G4int verbose)
{
  fVerbose = verbose;
}

inline G4double G4VMoleculeCounterInternalBase::GetActiveLowerBound() const
{
  return fActiveLowerBound;
}
inline void G4VMoleculeCounterInternalBase::SetActiveLowerBound(G4double time, G4bool inclusive)
{
  fActiveLowerBound = time;
  fActiveLowerBoundInclusive = inclusive;
}

inline G4double G4VMoleculeCounterInternalBase::GetActiveUpperBound() const
{
  return fActiveUpperBound;
}
inline void G4VMoleculeCounterInternalBase::SetActiveUpperBound(G4double time, G4bool inclusive)
{
  fActiveUpperBound = time;
  fActiveUpperBoundInclusive = inclusive;
}

inline G4bool G4VMoleculeCounterInternalBase::GetActiveLowerBoundInclusive() const
{
  return fActiveLowerBoundInclusive;
}
inline G4bool G4VMoleculeCounterInternalBase::GetActiveUpperBoundInclusive() const
{
  return fActiveUpperBoundInclusive;
}

inline G4bool G4VMoleculeCounterInternalBase::IsTimeBelowLowerBound(G4double time) const
{
  return (fActiveLowerBoundInclusive && time < fActiveLowerBound)
         || (!fActiveLowerBoundInclusive && time <= fActiveLowerBound);
}
inline G4bool G4VMoleculeCounterInternalBase::IsTimeAboveUpperBound(G4double time) const
{
  return (fActiveUpperBoundInclusive && time > fActiveUpperBound)
         || (!fActiveUpperBoundInclusive && time >= fActiveUpperBound);
}
inline G4bool G4VMoleculeCounterInternalBase::IsActiveAtGlobalTime(G4double time) const
{
  return !(IsTimeBelowLowerBound(time) || IsTimeAboveUpperBound(time));
}

inline G4bool G4VMoleculeCounterInternalBase::GetCheckTimeConsistencyWithScheduler() const
{
  return fCheckTimeIsConsistentWithScheduler;
}
inline void G4VMoleculeCounterInternalBase::SetCheckTimeConsistencyWithScheduler(G4bool flag)
{
  fCheckTimeIsConsistentWithScheduler = flag;
}

inline G4bool G4VMoleculeCounterInternalBase::GetCheckRecordedTimeConsistency() const
{
  return fCheckRecordedTimesAreConsistent;
}
inline void G4VMoleculeCounterInternalBase::SetCheckRecordedTimeConsistency(G4bool flag)
{
  fCheckRecordedTimesAreConsistent = flag;
}

inline const G4MoleculeCounterTimeComparer& G4VMoleculeCounterInternalBase::GetTimeComparer() const
{
  return fTimeComparer;
}

inline void G4VMoleculeCounterInternalBase::SetTimeComparer(const G4MoleculeCounterTimeComparer& comparer)
{
  if (fIsInitialized) {
    G4Exception("G4VMoleculeCounterInternalBase::SetTimeComparer()", "AlreadyInitialized", JustWarning,
                "Molecule counter was already initialized, assigning the time comparer now may "
                "have no effect!");
  }
  fTimeComparer = comparer;
}

//------------------------------------------------------------------------------

#endif
