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
// $Id: G3MatTableEntry.hh 67982 2013-03-13 10:36:03Z gcosmo $
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
    G3MatTableEntry& operator=(const G3MatTableEntry& right);
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
