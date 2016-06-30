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
//
// $Id: G4IsotopeProperty.hh 96314 2016-04-06 07:21:51Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Add Magnetic Moment                      14 Mar  07 H.Kurashige
//      Add isomer level              30 Apr. H.Kurashige


#ifndef G4IsotopeProperty_h
#define G4IsotopeProperty_h 1

#include "globals.hh"
class G4DecayTable;
#include "G4Ions.hh"
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

  // Set/Get Magentic Moment
  G4double      GetMagneticMoment() const;
  void          SetMagneticMoment(G4double M);

  // Set/Get Excited Energy
  G4double      GetEnergy() const;
  void          SetEnergy(G4double  E);

  // Set/Get isomer lervel
  G4int         GetIsomerLevel() const;
  void          SetIsomerLevel(G4int  level);

  // Set/Get floating level base
  G4Ions::G4FloatLevelBase GetFloatLevelBase() const;
  void                     SetFloatLevelBase(G4Ions::G4FloatLevelBase flb);
  void                     SetFloatLevelBase(G4int flbIndex);

  // Set/Get life time
  G4double      GetLifeTime() const;
  void          SetLifeTime(G4double  T);

  // Set/Get decay table
  G4DecayTable* GetDecayTable() const;
  void          SetDecayTable(G4DecayTable*  table);

  // Dump out information
  void          DumpInfo() const;

 private:
  G4int         fAtomicNumber;     // number of proton
  G4int         fAtomicMass;       // number of nucleon 
  G4int         fISpin;            // total angular momentum (in unit of 1/2)
  G4double      fEnergy;           // excited energy
  G4double      fLifeTime;         // lifeTime 
  G4DecayTable* fDecayTable;       // decay Table
  G4double      fMagneticMoment;   // magnetic moment 
  G4int         fIsomerLevel;      // isomer level 
  G4Ions::G4FloatLevelBase fFloatLevelBase; // floating level base
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
  G4double  G4IsotopeProperty::GetMagneticMoment() const
{
  return fMagneticMoment;
}

inline
  void     G4IsotopeProperty::SetMagneticMoment(G4double M)
{
  fMagneticMoment = M;
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
 G4int  G4IsotopeProperty::GetIsomerLevel() const
{
  return fIsomerLevel; 
}
inline  
 void   G4IsotopeProperty::SetIsomerLevel(G4int  level)
{
  fIsomerLevel = level;
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

inline 
 G4Ions::G4FloatLevelBase G4IsotopeProperty::GetFloatLevelBase() const
{ 
  return fFloatLevelBase; 
}

inline 
 void G4IsotopeProperty::SetFloatLevelBase(G4Ions::G4FloatLevelBase flb)
{
  fFloatLevelBase = flb;
}

inline 
 void G4IsotopeProperty::SetFloatLevelBase(G4int flbIndex)
{
  fFloatLevelBase = G4Ions::FloatLevelBase(flbIndex);
}

#endif








