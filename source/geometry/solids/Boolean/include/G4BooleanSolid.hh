// Abstract base class for solids created by boolean operations
// between other solids
//
// History:
//
// 10.09.98 V.Grichine, creation according J. Apostolakis's recommendations

#ifndef G4BOOLEANSOLID_HH
#define G4BOOLEANSOLID_HH

#include "G4VSolid.hh"
#include "G4DisplacedSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"


class G4BooleanSolid: public G4VSolid
{
public:
                  G4BooleanSolid( const G4String& pName,
                                        G4VSolid* pSolidA ,
                                        G4VSolid* pSolidB   ) ;

                  G4BooleanSolid( const G4String& pName,
                                        G4VSolid* pSolidA ,
                                        G4VSolid* pSolidB,
                                        G4RotationMatrix* rotMatrix,
                                  const G4ThreeVector& transVector    ) ;

                  G4BooleanSolid( const G4String& pName,
                                        G4VSolid* pSolidA ,
                                        G4VSolid* pSolidB , 
                                  const G4Transform3D& transform   ) ;

                 virtual ~G4BooleanSolid() ;


protected:
                  G4VSolid* fPtrSolidA ;
                  G4VSolid* fPtrSolidB ;

private:
                  G4bool  createdDisplacedSolid; // If & only if this object 
                                                // created it, it must delete it
} ;

#endif
