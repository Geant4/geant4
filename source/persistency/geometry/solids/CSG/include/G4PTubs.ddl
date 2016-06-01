// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PTubs.ddl,v 2.0 1998/07/02 16:13:11 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// class G4PTubs
//
// History:
// 19.06.98 A.Kimura Converted G4Tubs.hh

#ifndef G4PTUBS_DDL
#define G4PTUBS_DDL

#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Tubs;

class G4PTubs : public G4PCSGSolid {
public:
    G4PTubs(const G4Tubs* theTubs);
    virtual ~G4PTubs();

    G4VSolid* MakeTransientObject() const;

    // Naming method (pseudo-RTTI : run-time type identification)
    virtual G4GeometryType  GetEntityType() const { return G4String("G4Tubs"); }

protected:

    G4double fRMin,fRMax,fDz,fSPhi,fDPhi;

};
   	
#endif
