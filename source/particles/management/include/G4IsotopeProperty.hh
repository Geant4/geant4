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
// G4IsotopeProperty
//
// Class description:
//
// G4IsotopeProperty contains properties of an isotope

// Author: H.Kurashige, 5 October 1999
// --------------------------------------------------------------------
#ifndef G4IsotopeProperty_hh
#define G4IsotopeProperty_hh 1

#include "globals.hh"
#include "G4Ions.hh"

class G4DecayTable;
class G4IsotopeProperty
{
  public:

    G4IsotopeProperty();
    virtual ~G4IsotopeProperty();
      // Constructor & destructor

    G4IsotopeProperty(const G4IsotopeProperty& right);
      // Copy constructor  

    G4IsotopeProperty& operator=(G4IsotopeProperty& right);
      // Assignment operator
 
    G4bool operator==(const G4IsotopeProperty &right) const;
    G4bool operator!=(const G4IsotopeProperty &right) const;
      // Equality operators

    inline G4int GetAtomicNumber() const;
    inline void  SetAtomicNumber(G4int Z);
      // Set/Get Atomic Number

    inline G4int GetAtomicMass() const;
    inline void  SetAtomicMass(G4int A);
      // Set/Get Atomic Mass

    inline G4int GetiSpin() const;
    inline void  SetiSpin(G4int J);
      // Set/Get spin

    inline G4double GetMagneticMoment() const;
    inline void     SetMagneticMoment(G4double M);
      // Set/Get Magnetic Moment

    inline G4double GetEnergy() const;
    inline void     SetEnergy(G4double  E);
      // Set/Get Excited Energy

    inline G4int GetIsomerLevel() const;
    inline void  SetIsomerLevel(G4int level);
      // Set/Get isomer level

    inline G4Ions::G4FloatLevelBase GetFloatLevelBase() const;
    inline void SetFloatLevelBase(G4Ions::G4FloatLevelBase flb);
    inline void SetFloatLevelBase(G4int flbIndex);
      // Set/Get floating level base

    inline G4double GetLifeTime() const;
    inline void     SetLifeTime(G4double T);
      // Set/Get life time

    inline G4DecayTable* GetDecayTable() const;
    inline void SetDecayTable(G4DecayTable* table);
      // Set/Get decay table

    void DumpInfo() const;
      // Dump out information

  private:

    G4int         fAtomicNumber = 0; // number of proton
    G4int         fAtomicMass = 0;   // number of nucleon 
    G4int         fISpin = 0;        // total angular momentum (in unit of 1/2)
    G4double      fEnergy = 0.0;     // excited energy
    G4double      fLifeTime = -1.0;  // lifeTime 
    G4DecayTable* fDecayTable = nullptr;      // decay Table
    G4double      fMagneticMoment = 0.0;      // magnetic moment 
    G4int         fIsomerLevel = -1;          // isomer level 
    G4Ions::G4FloatLevelBase fFloatLevelBase; // floating level base
};

// ------------------------
// Inline methods
// ------------------------

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
void G4IsotopeProperty::SetIsomerLevel(G4int level)
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
