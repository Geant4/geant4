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
// $Id: G3MedTableEntry.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// The class associates the G3 tracking medium index
// with the corresponding G4Material, G4MagneticField, G4UserLimits 
// instances and sensitivity flag (isvol).

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3MEDTABLEENTRYH_HH
#define G3MEDTABLEENTRYH_HH 1

#include "globals.hh"

class G4Material;
class G4MagneticField;
class G4UserLimits;

class G3MedTableEntry 
{
  public:  // with description

    G3MedTableEntry(G4int id, G4Material* material, G4MagneticField* field,
       G4UserLimits* limits, G4int isvol);
    G3MedTableEntry(const G3MedTableEntry& right);
    virtual ~G3MedTableEntry();
    
    // operators
    G3MedTableEntry& operator=(const G3MedTableEntry& right);
    G4int operator==(const G3MedTableEntry& right) const;
    G4int operator!=(const G3MedTableEntry& right) const;

    // set methods
    void SetMaterial(G4Material* material);
    void SetField(G4MagneticField* field);
    void SetLimits(G4UserLimits* limits);
    void SetISVOL(G4int isvol);

    // get methods
    G4int GetID() const;
    G4Material* GetMaterial() const;
    G4MagneticField* GetField() const;
    G4UserLimits* GetLimits() const;
    G4int GetISVOL() const;
    
  private:

    // data members  
    G4int             fID;
    G4Material*       fMaterial;
    G4MagneticField*  fField;
    G4UserLimits*     fLimits;
    G4int             fISVOL;
    //G4double deemax;
    //G4double epsil;
};

// inline methods

inline void G3MedTableEntry::SetMaterial(G4Material* material)
{ fMaterial = material; }

inline void G3MedTableEntry::SetField(G4MagneticField* field)
{ fField = field; }

inline void G3MedTableEntry::SetLimits(G4UserLimits* limits)
{ fLimits = limits; }

inline void G3MedTableEntry::SetISVOL(G4int isvol)
{ fISVOL = isvol; }

inline G4int G3MedTableEntry::GetID() const
{ return fID; }

inline G4Material* G3MedTableEntry::GetMaterial() const
{ return fMaterial; }

inline G4MagneticField* G3MedTableEntry::GetField() const
{ return fField; }

inline G4UserLimits* G3MedTableEntry::GetLimits() const
{ return fLimits; }

inline G4int G3MedTableEntry::GetISVOL() const
{ return fISVOL; }

#endif
