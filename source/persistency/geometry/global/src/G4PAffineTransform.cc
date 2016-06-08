// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAffineTransform.cc,v 1.1.2.1 1999/12/07 20:50:06 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// class G4PAffineTransform
//
// History:
// 10.11.99 Y.Morita  Initial version

#include "G4PAffineTransform.hh"
#include "G4AffineTransform.hh"

G4PAffineTransform::G4PAffineTransform( const G4AffineTransform aTrans )
: rxx(aTrans[0]), rxy(aTrans[1]), rxz(aTrans[2]),
  ryx(aTrans[4]), ryy(aTrans[5]), ryz(aTrans[6]),
  rzx(aTrans[8]), rzy(aTrans[9]), rzz(aTrans[10]),
   tx(aTrans[12]), ty(aTrans[13]), tz(aTrans[14]) 
{;}

G4PAffineTransform::~G4PAffineTransform()
{;}

G4AffineTransform G4PAffineTransform::MakeTransientObject()
{
  G4RotationMatrix aRot;
  G4AffineTransform* aTrans = new G4AffineTransform
           ( aRot.rotateAxes(G4ThreeVector(rxx, ryx, rzx),
                             G4ThreeVector(rxy, ryy, rzy),
                             G4ThreeVector(rxz, ryz, rzz)),
                             G4ThreeVector( tx,  ty,  tz)  );
  G4AffineTransform trans = *aTrans;
  delete aTrans;
  return trans;
}

