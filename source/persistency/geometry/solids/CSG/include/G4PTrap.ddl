// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PTrap.ddl,v 1.1 1999/01/07 16:10:51 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PTrap
//
// History:
// 19.06.98 A.Kimura Converted G4Trap.hh

#ifndef G4PTrap_DDL
#define G4PTrap_DDL

#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Trap;

#include "G4ThreeVector.hh"
#include <rw/tvordvec.h>
typedef RWTValOrderedVector<G4ThreeVector> G4ThreeVectorList;

struct PTrapSidePlane
{
    G4double a,b,c,d;		// Normal unit vector (a,b,c)  and offset (d)
				// => Ax+By+Cz+D=0  
};

class G4PTrap : public G4PCSGSolid {

public:

    // The most general constructor for G4Trap     

    G4PTrap(const G4Trap* theTrap);
    virtual ~G4PTrap() ;

    // Access functions
    G4VSolid* MakeTransientObject() const;

    void SetAllParameters ( G4double pDz,
			    G4double pTheta,
			    G4double pPhi,
			    G4double pDy1,
			    G4double pDx1,
			    G4double pDx2,
			    G4double pAlp1,
			    G4double pDy2,
			    G4double pDx3,
			    G4double pDx4,
			    G4double pAlp2);
	                              
    // Naming method (pseudo-RTTI : run-time type identification
    virtual G4GeometryType  GetEntityType() const { return G4String("G4Trap"); }

protected:

    G4bool MakePlanes();
    G4bool MakePlane( const G4ThreeVector& p1,
                      const G4ThreeVector& p2,
		      const G4ThreeVector& p3, 
		      const G4ThreeVector& p4,
		      PTrapSidePlane& plane ) ;

private:

    G4double fDz,fTthetaCphi,fTthetaSphi;
    G4double fDy1,fDx1,fDx2,fTalpha1;
    G4double fDy2,fDx3,fDx4,fTalpha2;
    PTrapSidePlane fPlanes[4];
};

#endif


//  **************************** End of G4Trap.hh *****************************************
