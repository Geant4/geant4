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
#ifndef G4NuclideTable_h
#define G4NuclideTable_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4NuclideTable.hh
//
// Date:                10/10/13
// Author:              T.Koi
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HISTORY
// Based on G4IsomerTable
//
///////////////////////////////////////////////////////////////////////////////

//
#include "globals.hh"
#include "G4IsotopeProperty.hh"
#include "G4VIsotopeTable.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"

#include <vector>
#include <cmath>

class G4NuclideTableMessenger;

////////////////////////////////////////////////////////////////////////////////
//
class G4NuclideTable : public G4VIsotopeTable
{
  // class description
  //   G4NuclideTable is the table of pointers to G4IsotopeProperty,
  //   which has magnetic moment and spin.
  //   Data File name is given by G4ENSDFSTATEDATA
  //   
  //
   private:
      G4NuclideTable();
      //G4NuclideTable ( const G4NuclideTable& p ){};
      //G4NuclideTable& operator=( const G4NuclideTable& p ){};

   public:
      static G4NuclideTable* GetInstance(); 
      static G4NuclideTable* GetNuclideTable() { return GetInstance(); };
public:
  typedef std::vector<G4IsotopeProperty*> G4IsotopeList;

protected:
  void FillHardCodeList();

public:
  // destructor
  virtual ~G4NuclideTable();

public:

  //
  void GenerateNuclide();

  void SetThresholdOfHalfLife( G4double );
  G4double GetThresholdOfHalfLife() { return threshold_of_half_life; };

  void SetLevelTolerance( G4double x ) { flevelTolerance=x;};
  G4double GetLevelTolerance() { return flevelTolerance; };

  void AddState(G4int,G4int,G4double,G4double,G4int ionJ=0,G4double ionMu=0.0);
  void AddState(G4int,G4int,G4double,G4int,G4double,G4int ionJ=0,G4double ionMu=0.0);
  void AddState(G4int,G4int,G4double,G4Ions::G4FloatLevelBase,G4double,G4int ionJ=0,G4double ionMu=0.0);
   
  //G4bool Exists(G4double,G4double,G4double);
  
  //Use GetIsotope(G4int Z, G4int A, G4double E)
  //void FillProperty(G4ParticleDefinition*);

  size_t GetSizeOfIsotopeList()
  { return ( fIsotopeList ? fIsotopeList->size() : static_cast<size_t>(0)) ; };
  
  //
  // with description
  //
  virtual G4IsotopeProperty* GetIsotope(G4int Z, G4int A, G4double E,
            G4Ions::G4FloatLevelBase flb=G4Ions::G4FloatLevelBase::no_Float);
  virtual G4IsotopeProperty* GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl=0);
  //
  //   again it will replace the pure virtual one in the abstract base class.
  //
  //   Z: Atomic Number
  //   A: Atomic Mass
  //   E: Excitaion energy
  //   flb: floating level base (enum defined in G4Ions.hh)
  //    or
  //   lvl: isomer level
  //  

  size_t  entries() const; 
  G4IsotopeProperty* GetIsotopeByIndex(size_t idx) const;

public:
  // utility methods
  static inline G4double GetTrancationError( G4double eex )
  { G4double tolerance= G4NuclideTable::GetInstance()->GetLevelTolerance();
  return eex - (G4long)(eex/tolerance)*tolerance; }
  static inline G4double Round( G4double eex )
  { G4double tolerance= G4NuclideTable::GetInstance()->GetLevelTolerance();
    return round(eex/tolerance)*tolerance; }
  static inline G4long Trancate( G4double eex )
  { G4double tolerance= G4NuclideTable::GetInstance()->GetLevelTolerance();
    return (G4long)(eex/tolerance); }
  static inline G4double Tolerance()
  { return G4NuclideTable::GetInstance()->GetLevelTolerance(); }

private:

  G4double threshold_of_half_life;         //threshold values of half-life of current run
  G4double minimum_threshold_of_half_life; //The minimum value of threshold values of half-life during entire runs
  G4IsotopeList*        fUserDefinedList;

  //Design Change on v10.02
  //pre_load_list: Having state data for current run defined by threshold_of_half_life
  //full_list:Keeping all state data during  running application defined by minimum_threshold_of_half_life
  //       ionCode                Ex. Energy
  std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > > map_pre_load_list;
  std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > > map_full_list;

  // Table of Nuclide Property
  //  0: Z
  //  1: A 
  //  2: Energy [keV]
  //  3: Life Time [ns]
  //  4: Spin  [h_bar/2]
  //  5: Magnetic Moment [joule/tesla]
  enum {idxZ=0, idxA,idxEnergy, idxLife, idxSpin, idxMu };

  G4IsotopeList*        fIsotopeList;
  G4double flevelTolerance;
  G4NuclideTableMessenger* fMessenger;

private:
  G4double StripFloatLevelBase(G4double E, G4int& flbIndex);
  G4Ions::G4FloatLevelBase StripFloatLevelBase( G4String );
};

inline
 size_t  G4NuclideTable::entries() const
{
  return (fIsotopeList ? fIsotopeList->size() : static_cast<size_t>(0) );
}

inline
  G4IsotopeProperty* G4NuclideTable::GetIsotopeByIndex(size_t idx) const
{
  if ( fIsotopeList && idx<fIsotopeList->size()) return (*fIsotopeList)[idx];
  else                          return 0;
}

#endif
