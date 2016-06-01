// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PBox.ddl,v 2.0 1998/07/02 16:12:56 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4PBox
//
// History:
// 19.06.98 A.Kimura Converted from G4Box.hh

#ifndef G4PBOX_DDL
#define G4PBOX_DDL

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

