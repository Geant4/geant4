// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polymarker.hh,v 1.2 1999-05-19 08:33:42 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  November 1996

#ifndef G4POLYMARKER_HH
#define G4POLYMARKER_HH

#include "G4VMarker.hh"
#include "G4Point3DList.hh"

class G4Polymarker: public G4VMarker, public G4Point3DList {
public:
  friend ostream& operator << (ostream& os, const G4Polymarker& marker);
  enum MarkerType {line, dots, circles, squares};
  G4Polymarker ();
  virtual G4Visible&  operator = (const G4Visible& right);
  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
  virtual G4VMarker&  operator = (const G4VMarker& right);
  MarkerType GetMarkerType () const;
  void SetMarkerType (MarkerType type);
private:
  MarkerType fMarkerType;
};

#include "G4Polymarker.icc"

#endif
