// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PCons.ddl,v 2.0 1998/07/02 16:12:58 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4Cons
//
// History:
// 19.06.98 A.Kimura Converted G4Cons.hh

#ifndef G4PCons_DDL
#define G4PCons_DDL

#include "G4PCSGSolid.hh"

class G4Cons;
class G4VSolid;

class G4PCons : public G4PCSGSolid {
public:
    G4PCons(const G4Cons* theCons);
    virtual ~G4PCons() ;

    G4VSolid* MakeTransientObject() const;

    // Naming method (pseudo-RTTI : run-time type identification
    virtual G4GeometryType  GetEntityType() const {return G4String("G4Cons");}

private:

    G4double fRmin1, fRmin2,
             fRmax1, fRmax2,
             fDz,
	     fSPhi, fDPhi;
};

#endif
