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
#ifndef G4IsotopeMagneticMomentTable_h
#define G4IsotopeMagneticMomentTable_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4IsotopeMagneticMomentTable.hh
//
// Date:                16/03/07
// Author:              H.Kurashige
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HISTORY
////////////////////////////////////////////////////////////////////////////////// IsomerLevel is added                         30 Apr. 2013  H.Kurashige

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
class G4IsotopeMagneticMomentTable : public G4VIsotopeTable
{
  // class description
  //   G4IsotopeMagneticMomentTable is the table of pointers to G4IsotopeProperty,
  //   which has magnetic moment and spin.
  //   Data File name is given by G4IONMAGNETICMOMENT 
  //   
public:
  //
  typedef std::vector<G4IsotopeProperty*> G4IsotopeList;
  typedef std::vector<G4String>           G4IsotopeNameList;

public:
  // constructor
  //
  G4IsotopeMagneticMomentTable();

protected:
  // hide copy construictor and assignment operator as protected
  G4IsotopeMagneticMomentTable(const  G4IsotopeMagneticMomentTable &right);
  G4IsotopeMagneticMomentTable & operator= (const  G4IsotopeMagneticMomentTable &right);

public:
  // destructor
  virtual ~G4IsotopeMagneticMomentTable();

public:
  // with description
  //
  virtual G4bool FindIsotope(G4IsotopeProperty* property);
  // The FindIsotope function will replace the pure virtual one defined in the
  // abstract base class G4VIstopeTable.  We don't use this fuction in this
  // implementation, instead we use the next function.
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
  //          -- currently ignored
  //    or
  //    G4int  level: isomer level
  // 

private:
  G4int GetVerboseLevel() const;
  // get Verbose Level defined in G4ParticleTable

private:

  G4IsotopeList         fIsotopeList;

  static const G4double levelTolerance;
  static const G4double nuclearMagneton;
};


inline 
 G4int G4IsotopeMagneticMomentTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

#endif






