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
// $Id: G4PBox.ddl,v 1.3.4.1 2001/06/28 19:11:31 gunter Exp $
// GEANT4 tag $Name:  $
//
// class G4PBox
//
// History:
// 19.06.98 A.Kimura Converted from G4Box.hh

#ifndef G4PBOX_DDL
#define G4PBOX_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Box;

class G4PBox : public G4PCSGSolid {
public:
    G4PBox(const G4Box* theBox);
    virtual ~G4PBox();

    G4VSolid* MakeTransientObject() const;

    // Naming method (pseudo-RTTI : run-time type identification
    virtual G4GeometryType  GetEntityType() const { return G4String("G4Box"); }

private:
    G4double fDx,fDy,fDz;
};

#endif

