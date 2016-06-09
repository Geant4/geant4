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
#ifndef G4RIsotopeTable_h
#define G4RIsotopeTable_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4RIsotopeTable.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4IsotopeProperty.hh"
#include "G4VIsotopeTable.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"

#include <vector>

////////////////////////////////////////////////////////////////////////////////
//
class G4RIsotopeTable : public G4VIsotopeTable
{
  // class description
  //   G4RIsotopeTable is the table of pointers to G4IsotopeProperty.
  //   It is derived from the abstract base class G4VIsotopeTable
  // class description - end
public:
  //
  typedef std::vector<G4IsotopeProperty*> G4IsotopeList;
  typedef std::vector<G4String> G4IsotopeNameList;

public:
  // constructor
  G4RIsotopeTable();

protected:
  // Hide copy constructor and assignment operator as protected
  G4RIsotopeTable(const G4RIsotopeTable& right);
  G4RIsotopeTable& operator = (const G4RIsotopeTable& right);

public:
  // destructor
  virtual ~G4RIsotopeTable();

public:
  // with description
  //
  virtual G4bool FindIsotope(G4IsotopeProperty* property);
  // The FindIsotope function will replace the pure virtual one defined in the
  // abstract base class G4VIstopeTable.  We don't use this fuction in this
  // implementation, instead we use the next function.
  //
  virtual G4IsotopeProperty* GetIsotope(G4int Z, G4int A, G4double E);
  //
  //   again it will replace the pure virtual one in the abstract base class.
  //
  //   Z: Atomic Number
  //   A: Atomic Mass
  //   E: Excitaion energy


  void AddUserDecayDataFile(G4int Z, G4int A,G4String filename);
  //Allow the user to replace the radio-active decay data provided in Geant4
  // by its own data file for a given isotope


private:

  G4String GetIsotopeName(G4int Z, G4int A, G4double E);
  //
  //   To generate a name for the isotope
  G4double GetMeanLifeTime (G4int Z, G4int A, G4double& E);
  // to get the mean life time and its decaytable from the database
  //
  G4IsotopeProperty* GetIsotope(G4int index) const;
  // Return the pointer of index-th isotope in the table

  G4int GetVerboseLevel() const;
  // get Verbose Level defined in G4ParticleTable

  G4int Entries() const;

private:

  G4IsotopeList         fIsotopeList;
  G4IsotopeNameList     fIsotopeNameList;
  static const G4double levelTolerance;

  //User define radioactive decay data files replacing some files in the G4RADECAY database
  std::map<G4int, G4String> theUserRadioactiveDataFiles;

};

inline G4int G4RIsotopeTable::Entries() const
  {return fIsotopeNameList.size();}

inline G4IsotopeProperty*  G4RIsotopeTable::GetIsotope(G4int index) const
{
  if ( (index >=0 ) && (index < Entries()) ){
    return fIsotopeList[index];
  } else {
    return 0; 
  } 
}

#endif






