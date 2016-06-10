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
// $Id: G3RotTableEntry.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// The class associates the G3 rotation matrix index to its
// corresponding G4RotationMatrix object.

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3ROTTABLEENTRY_HH
#define G3ROTTABLEENTRY_HH 1

#include "globals.hh"
#include "G4RotationMatrix.hh"

class G3toG4RotationMatrix;

class G3RotTableEntry 
{
  public:  // with description

    G3RotTableEntry(G4int id, G4RotationMatrix* matrix);
    G3RotTableEntry(const G3RotTableEntry& right);
    virtual ~G3RotTableEntry();
    
    // operators
    G3RotTableEntry& operator=(const G3RotTableEntry& right);
    G4int operator==(const G3RotTableEntry& right) const;
    G4int operator!=(const G3RotTableEntry& right) const;

    // get methods
    G4int       GetID() const;
    G4RotationMatrix* GetMatrix() const;
    
  private:

    // data members  
    G4int              fID;
    G4RotationMatrix*  fMatrix;
};

// inline methods

inline G4int G3RotTableEntry::GetID() const
{ return fID; }

inline G4RotationMatrix* G3RotTableEntry::GetMatrix() const
{ return fMatrix; }

#endif
