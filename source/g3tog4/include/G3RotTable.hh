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
// $Id: G3RotTable.hh,v 1.15 2006-06-29 18:12:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// G3 rotations table.
// Maps G3 rotations indices to their G4RotationMatrix object counterparts.

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3ROTTABLEH_HH
#define G3ROTTABLEH_HH 1

#include "G3RotTableEntry.hh"

#include "globals.hh"

#include <vector>

class G4Material;

typedef std::vector<G3RotTableEntry*>  G3RotMatrixVector;

class G3RotTable
{
  public:  // with description

    G3RotTable();
    virtual ~G3RotTable();
    
    // methods

    G4RotationMatrix* Get(G4int id) const;
    void Put(G4int id, G4RotationMatrix* matrix);
    void Clear();

  private:

    G3RotMatrixVector*  fRotVector;
};

extern G3RotTable G3Rot;

#endif
