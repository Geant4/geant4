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
// $Id: G3MatTableEntry.hh,v 1.5 2001-07-11 09:58:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// The class associates the G3 material index with the
// corresponding G4Material object.

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3MATTABLEENTRY_HH
#define G3MATTABLEENTRY_HH 1

#include "globals.hh"

class G4Material;

class G3MatTableEntry 
{
  public: // with description
  
    G3MatTableEntry(G4int id, G4Material* material);
    G3MatTableEntry(const G3MatTableEntry& right);
    virtual ~G3MatTableEntry();
    
    // operators
    const G3MatTableEntry& operator=(const G3MatTableEntry& right);
    G4int operator==(const G3MatTableEntry& right) const;
    G4int operator!=(const G3MatTableEntry& right) const;

    // get methods
    G4int       GetID() const;
    G4Material* GetMaterial() const;
    
  private:
  
    // data members  
    G4int        fID;
    G4Material*  fMaterial;
};

// inline methods

inline G4int G3MatTableEntry::GetID() const
{ return fID; }

inline G4Material* G3MatTableEntry::GetMaterial() const
{ return fMaterial; }

#endif
