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
// G4NuclideTable
//
// Class description:
//
// Table of pointers to G4IsotopeProperty, which has magnetic moment
// and spin. Data File name is given by G4ENSDFSTATEDATA.
// Table based on G4IsomerTable.

// Author: T.Koi, SLAC - 10 October 2013
// --------------------------------------------------------------------
#ifndef G4NuclideTable_hh
#define G4NuclideTable_hh 1

#include <vector>
#include <cmath>

#include "globals.hh"
#include "G4IsotopeProperty.hh"
#include "G4VIsotopeTable.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"

class G4NuclideTableMessenger;

class G4NuclideTable : public G4VIsotopeTable
{
  public:

    using G4IsotopeList = std::vector<G4IsotopeProperty*>;

    virtual ~G4NuclideTable();

    G4NuclideTable (const G4NuclideTable&) = delete;
    G4NuclideTable& operator = (const G4NuclideTable&) = delete;

    static G4NuclideTable* GetInstance(); 
    static G4NuclideTable* GetNuclideTable();

    void GenerateNuclide();

    void SetThresholdOfHalfLife(G4double);
    inline G4double GetThresholdOfHalfLife();

    void SetMeanLifeThreshold(G4double);
    inline G4double GetMeanLifeThreshold();

    inline void SetLevelTolerance(G4double x);
    inline G4double GetLevelTolerance();

    void AddState(G4int, G4int, G4double, G4double,
                  G4int ionJ=0, G4double ionMu=0.0);
    void AddState(G4int, G4int, G4double, G4int, G4double,
                  G4int ionJ=0, G4double ionMu=0.0);
    void AddState(G4int, G4int, G4double, G4Ions::G4FloatLevelBase, G4double,
                  G4int ionJ=0, G4double ionMu=0.0);

    inline std::size_t GetSizeOfIsotopeList();
  
    virtual G4IsotopeProperty* GetIsotope(G4int Z, G4int A, G4double E,
            G4Ions::G4FloatLevelBase flb=G4Ions::G4FloatLevelBase::no_Float);
    virtual G4IsotopeProperty* GetIsotopeByIsoLvl(G4int Z, G4int A,
            G4int lvl=0);
      // It will replace the pure virtual one in the abstract base class.
      //   Z: Atomic Number
      //   A: Atomic Mass
      //   E: Excitaion energy
      //   flb: floating level base (enum defined in G4Ions.hh)
      //    or
      //   lvl: isomer level

    inline std::size_t entries() const; 
    inline G4IsotopeProperty* GetIsotopeByIndex(std::size_t idx) const;

    // utility methods

    static G4double GetTruncationError(G4double eex);
    static G4double Round(G4double eex);
    static G4long Truncate(G4double eex);
    static G4double Tolerance();

  private:

    G4NuclideTable();

    G4double StripFloatLevelBase(G4double E, G4int& flbIndex);
    G4Ions::G4FloatLevelBase StripFloatLevelBase(const G4String&);

  private:

    G4double mean_life_threshold = 0.0;
    G4double minimum_mean_life_threshold = DBL_MAX;

    G4IsotopeList* fUserDefinedList = nullptr;

    std::map<G4int, std::multimap<G4double, G4IsotopeProperty*> > map_pre_load_list;
      // pre_load_list: contains state data for current run defined
      // by mean_life_threshold
    std::map<G4int, std::multimap<G4double, G4IsotopeProperty*> > map_full_list;
      // full_list: keeps all state data during running application
      // defined by minimum_mean_life_threshold

    enum { idxZ=0, idxA, idxEnergy, idxLife, idxSpin, idxMu };
      // Table of Nuclide Property
      //  0: Z
      //  1: A 
      //  2: Energy [keV]
      //  3: Life Time [ns]
      //  4: Spin  [h_bar/2]
      //  5: Magnetic Moment [joule/tesla]

    G4IsotopeList* fIsotopeList = nullptr;
    G4double flevelTolerance = 0.0;
    G4NuclideTableMessenger* fMessenger = nullptr;
};

// ------------------------
// Inline methods
// ------------------------

inline
G4double G4NuclideTable::GetThresholdOfHalfLife()
{
  return mean_life_threshold*0.69314718;
}

inline
G4double G4NuclideTable::GetMeanLifeThreshold()
{
  return mean_life_threshold;
}

inline
void G4NuclideTable::SetLevelTolerance(G4double x)
{
  flevelTolerance = x;
}

inline
G4double G4NuclideTable::GetLevelTolerance()
{
  return flevelTolerance;
}

inline
std::size_t G4NuclideTable::GetSizeOfIsotopeList()
{
  return (fIsotopeList ? fIsotopeList->size() : static_cast<size_t>(0) );
}

inline
std::size_t G4NuclideTable::entries() const
{
  return (fIsotopeList ? fIsotopeList->size() : std::size_t(0) );
}

inline
G4IsotopeProperty* G4NuclideTable::GetIsotopeByIndex(std::size_t idx) const
{
  if (fIsotopeList && idx<fIsotopeList->size() ) return (*fIsotopeList)[idx];
  else return nullptr;
}

#endif
