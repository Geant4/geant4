// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FieldTrack.cc,v 1.2 1999-12-15 14:49:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4FieldTrack.hh"

G4std::ostream& operator<<( G4std::ostream& os, G4FieldTrack& SixVec)
{
     G4double *SixV = SixVec.SixVector;
     os << " X= " << SixV[0] << " " << SixV[1] << " " << SixV[2] << " ";
     os << " V= " << SixV[3] << " " << SixV[4] << " " << SixV[5] << " ";
     os << " l= " << SixVec.CurveS();
     return os;
}
