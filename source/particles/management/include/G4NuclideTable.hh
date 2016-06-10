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

  void SetThresholdOfHalfLife( G4double t ) { threshold_of_half_life=t;};
  G4double GetThresholdOfHalfLife() { return threshold_of_half_life; };

  void AddState(G4int,G4int,G4double,G4double,G4int,G4double);

   
  //G4bool Exists(G4double,G4double,G4double);
  
  //Use GetIsotope(G4int Z, G4int A, G4double E)
  //void FillProperty(G4ParticleDefinition*);

  size_t GetSizeOfIsotopeList(){ return ( fIsotopeList ? fIsotopeList->size() : static_cast<size_t>(0)) ; };
  
  //
  // with description
  //
  virtual G4IsotopeProperty* GetIsotope(G4int Z, G4int A, G4double E);
  virtual G4IsotopeProperty* GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl=0);
  //
  //   again it will replace the pure virtual one in the abstract base class.
  //
  //   Z: Atomic Number
  //   A: Atomic Mass
  //   E: Excitaion energy
  //    or
  //   lvl: isomer level
  //  

  size_t  entries() const; 
  G4IsotopeProperty* GetIsotopeByIndex(size_t idx) const;

public:
  // utility methods
  static inline G4double GetTrancationError( G4double eex )
  { return eex - (G4long)(eex/levelTolerance)*levelTolerance; }
  static inline G4double Round( G4double eex )
  { return (G4long)(eex/levelTolerance)*levelTolerance; }
  static inline G4long Trancate( G4double eex )
  { return (G4long)(eex/levelTolerance); }
  static inline G4double Tolerance()
  { return levelTolerance; }

private:

  G4double threshold_of_half_life;
  G4IsotopeList*        fUserDefinedList;

  std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > > map_pre_load_list;
  std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > > map_hard_code_list;
  std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > > map_full_list;

  //enum {nEntries=3075,MaxA=260, MinZ=2, MaxZ=100};
  enum {nEntries_ground_state=2910};
  enum {nEntries_excite_state=3898};
  //static const G4double isomerTable[nEntries][5];
  static const G4double groundStateTable[nEntries_ground_state][6];
  static const G4double exciteStateTable[nEntries_excite_state][6];
  // Table of Nuclide Property
  //  0: Z
  //  1: A 
  //  2: Energy [keV]
  //  3: Life Time [ns]
  //  4: Spin  [h_bar/2]
  //  5: Magnetic Moment [joule/tesla]
  enum {idxZ=0, idxA,idxEnergy, idxLife, idxSpin, idxMu };

  G4IsotopeList*        fIsotopeList;
  static const G4double levelTolerance;
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
