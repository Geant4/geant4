//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4IsotopeProperty.hh,v 1.3 2001-07-11 10:01:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige

#ifndef G4IsotopeProperty_h
#define G4IsotopeProperty_h 1

#include "globals.hh"
class    G4DecayTable;
class G4IsotopeProperty
{
 // Class Description
 //   G4IsotopeProperty contains properties of an isotope
 //

 public:
  G4IsotopeProperty();

  // copy construictor  
  G4IsotopeProperty(const  G4IsotopeProperty& right);

  // Assignment operator
  G4IsotopeProperty & operator=(G4IsotopeProperty& right);
 
  // equal / unequal operator
  G4int operator==(const G4IsotopeProperty &right) const;
  G4int operator!=(const G4IsotopeProperty &right) const;

  // destructor
  virtual ~G4IsotopeProperty();


 public:  // With Description
  // Set/Get Atomic Number
  G4int         GetAtomicNumber() const;
  void          SetAtomicNumber(G4int Z);

  // Set/Get Atomic Mass
  G4int         GetAtomicMass() const;
  void          SetAtomicMass(G4int A);

  // Set/Get spin
  G4int         GetiSpin() const;
  void          SetiSpin(G4int J);

  // Set/Get Excited Energy
  G4double      GetEnergy() const;
  void          SetEnergy(G4double  E);

  // Set/Get life time
  G4double      GetLifeTime() const;
  void          SetLifeTime(G4double  T);

  // Set/Get decay table
  G4DecayTable* GetDecayTable() const;
  void          SetDecayTable(G4DecayTable*  table);

  // Dump out information
  void          DumpInfo() const;

 private:
  G4int         fAtomicNumber; // number of proton
  G4int         fAtomicMass;   // number of nucleon 
  G4int         fISpin;        // total angular momentum (in unit of 1/2)
  G4double      fEnergy;       // excited energy
  G4double      fLifeTime;     // lifeTime 
  G4DecayTable* fDecayTable;   // decay Table
};

inline 
 G4int G4IsotopeProperty::GetAtomicNumber() const
{
  return fAtomicNumber;
}

inline 
 void G4IsotopeProperty::SetAtomicNumber(G4int Z)
{
    fAtomicNumber = Z;
}

inline 
 G4int G4IsotopeProperty::GetAtomicMass() const
{
  return fAtomicMass;
}

inline 
 void G4IsotopeProperty::SetAtomicMass(G4int A)
{
    fAtomicMass = A;
}

inline 
 G4int G4IsotopeProperty::GetiSpin() const
{
  return fISpin;
}

inline 
 void G4IsotopeProperty::SetiSpin(G4int J)
{
    fISpin = J;
}

inline 
 G4double G4IsotopeProperty::GetEnergy() const
{
  return fEnergy;
}

inline 
 void G4IsotopeProperty::SetEnergy(G4double E)
{
    fEnergy = E;
}

inline 
 G4double G4IsotopeProperty::GetLifeTime() const
{
  return fLifeTime;
}

inline 
 void G4IsotopeProperty::SetLifeTime(G4double T)
{
    fLifeTime = T;
}

inline 
 G4DecayTable* G4IsotopeProperty::GetDecayTable() const
{
  return fDecayTable;
}

inline 
 void G4IsotopeProperty::SetDecayTable(G4DecayTable* table)
{
    fDecayTable = table;
}
#endif










