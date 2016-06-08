// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PSphere.ddl,v 1.4 1999/12/15 14:51:25 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// class G4PSphere
//
// History:
// 19.06.98 A.Kimura Converted G4Sphere.hh

#ifndef G4PSphere_DDL
#define G4PSphere_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Sphere;

#include "G4ThreeVector.hh"
#include "g4rw/tvordvec.h"
typedef G4RWTValOrderedVector<G4ThreeVector> G4ThreeVectorList;

class G4PSphere : public G4PCSGSolid {
public:
    G4PSphere(const G4Sphere* theSphere);
		   
    virtual ~G4PSphere() ;

    G4VSolid* MakeTransientObject() const;

    // Naming method (pseudo-RTTI : run-time type identification
    virtual G4GeometryType  GetEntityType() const {return G4String("G4Sphere");}
       
private:

    G4double fRmin,fRmax,
             fSPhi,fDPhi,
	     fSTheta,fDTheta;
};
   	
#endif

