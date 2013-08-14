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
#ifndef G4IsomerTable_h
#define G4IsomerTable_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4IsomerTable.hh
//
// Date:                5/05/13
// Author:              H.Kurashige
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HISTORY
// 25/07/13   isomerTable is chaned to const member 
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
class G4IsomerTable : public G4VIsotopeTable
{
  // class description
  //   G4IsomerTable is the table of pointers to G4IsotopeProperty,
  //   which has magnetic moment and spin.
  //   Data File name is given by G4IONMAGNETICMOMENT 
  //   
public:
  //
  typedef std::vector<G4IsotopeProperty*> G4IsotopeList;

public:
  // constructor
  //
  G4IsomerTable();

protected:
  // hide copy construictor and assignment operator as protected
  G4IsomerTable(const  G4IsomerTable &right);
  G4IsomerTable & operator= (const  G4IsomerTable &right);

public:
  // destructor
  virtual ~G4IsomerTable();

public:
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

protected: 
  void FillIsotopeList();

private:
  enum {nEntries=3075,MaxA=260, MinZ=2, MaxZ=100};
  static const G4double isomerTable[nEntries][5];
  // Table of Isomer Property
  //  0: PID = Z*10000 + 10*A + Lvl
  //  1: Energy [keV]
  //  2: Life Time [ns]
  //  3: Spin  [h_bar/2]
  //  4: Magnetic Moment [joule/tesla]
  enum {idxPID=0, idxEnergy, idxLife, idxSpin, idxMu };

  G4IsotopeList*        fIsotopeList;
  static const G4double levelTolerance;
};


inline
 size_t  G4IsomerTable::entries() const
{
  return fIsotopeList->size();
}

inline
  G4IsotopeProperty* G4IsomerTable::GetIsotopeByIndex(size_t idx) const
{
  if (idx<fIsotopeList->size()) return (*fIsotopeList)[idx];
  else                          return 0;
}
#endif






