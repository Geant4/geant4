// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPara.ddl,v 1.1 1999/01/07 16:10:51 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PPara
//
// History:
// 19.06.98 A.Kimura Converted G4Para.hh


#ifndef G4PPara_DDL
#define G4PPara_DDL

#include "G4PCSGSolid.hh"

class G4Para;
class G4VSolid;

#include "G4ThreeVector.hh"
#include <rw/tvordvec.h>
typedef RWTValOrderedVector<G4ThreeVector> G4ThreeVectorList;

class G4PPara : public G4PCSGSolid {
public:

    G4PPara(const G4Para* thePara);
    virtual ~G4PPara() ;
    
				// Access functions
    G4VSolid* MakeTransientObject() const ;

				// Get  functions

    G4ThreeVector GetSymAxis() const {
	G4double cosTheta
	    = 1.0/sqrt(1+fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi) ;
     
	return G4ThreeVector(fTthetaCphi*cosTheta,fTthetaSphi*cosTheta,cosTheta) ;
   }

    // Naming method (pseudo-RTTI : run-time type identification

    virtual G4GeometryType  GetEntityType() const {return G4String("G4Para");}

private:
    G4double fDx,fDy,fDz;
    G4double fTalpha,fTthetaCphi,fTthetaSphi;
};
   	
#endif
