// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Text.hh,v 1.1 1999-01-07 16:09:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  17/11/96.

#ifndef G4TEXT_HH
#define G4TEXT_HH

#include "G4VMarker.hh"
#include "globals.hh"

class G4Text: public G4VMarker {
public:
  enum Layout {left, centre, right};
  G4Text (const G4String& text);
  G4Text (const G4String& text, const G4Point3D& pos);
  G4Text (const G4VMarker& marker);
  G4String GetText   () const;
  Layout   GetLayout () const;

  G4double GetXOffset () const ;
  G4double GetYOffset () const ;

  void SetText   (const G4String& text);
  void SetLayout (Layout layout);

  void SetOffset ( double dx,  double dy ) ;   

private:
  G4String fText;
  Layout   fLayout;
  G4double fXOffset, fYOffset ;
};

#include "G4Text.icc"

#endif
