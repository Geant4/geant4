// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polymarker.hh,v 1.5 1999-12-15 14:50:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  November 1996

// Class Description:
// A set of markers.
// Class Description - End:


#ifndef G4POLYMARKER_HH
#define G4POLYMARKER_HH

#include "G4VMarker.hh"
#include "G4Point3DList.hh"

class G4Polymarker: public G4VMarker, public G4Point3DList {

public: // With description

  friend G4std::ostream& operator << (G4std::ostream& os, const G4Polymarker& marker);
  enum MarkerType {line, dots, circles, squares};
  G4Polymarker ();
  virtual G4Visible&    operator = (const G4Visible& right);
  virtual G4VVisPrim&   operator = (const G4VVisPrim& right);
  virtual G4VMarker&    operator = (const G4VMarker& right);
  virtual G4Polymarker& operator = (const G4Polymarker& right);
  MarkerType GetMarkerType () const;
  void SetMarkerType (MarkerType type);
private:
  MarkerType fMarkerType;
};

#include "G4Polymarker.icc"

#endif
