// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisAttributes.hh,v 1.1 1999-01-07 16:09:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  23rd October 1996

#ifndef __G4VISATTRIBUTES_HH__
#define __G4VISATTRIBUTES_HH__

#include "G4ios.hh"
#include "globals.hh"
#include "G4Colour.hh"
#include "G4Color.hh"

class G4VisAttributes {

  friend ostream& operator << (ostream& os, const G4VisAttributes& a);
  friend G4bool   operator != (const G4VisAttributes& a1,
			       const G4VisAttributes& a2);

public:

  // Constructors, etc. - begin snippet.
  enum LineStyle {unbroken, dashed, dotted};
  enum ForcedDrawingStyle {wireframe, solid};

  G4VisAttributes ();
  G4VisAttributes (G4bool visibility);
  G4VisAttributes (const G4Colour& colour);
  G4VisAttributes (G4bool visibility, const G4Colour& colour);

  static const G4VisAttributes Invisible;
  // Constructors - end snippet.

  G4bool          IsVisible                      () const;
  G4bool          IsDaughtersInvisible           () const;
  const G4Colour& GetColour                      () const;
  const G4Color&  GetColor                       () const;
  LineStyle       GetLineStyle                   () const;
  G4double        GetLineWidth                   () const;
  G4bool          IsForceDrawingStyle            () const;
  ForcedDrawingStyle GetForcedDrawingStyle () const;

  // Set methods - begin snippet.
  void SetVisibility         (G4bool);
  void SetDaughtersInvisible (G4bool);
  void SetColour             (const G4Colour&);
  void SetColor              (const G4Color&);
  void SetColour             (G4double red, G4double green, G4double blue,
			      G4double alpha = 1.);
  void SetColor              (G4double red, G4double green, G4double blue,
			      G4double alpha = 1.);
  void SetLineStyle          (LineStyle);
  void SetLineWidth          (G4double);
  void SetForceWireframe     (G4bool);
  void SetForceSolid         (G4bool);
  // Set methods - end snippet.

private:

  // Available attributes - begin snippet.
  G4bool      fVisible;            // Visibility flag
  G4bool      fDaughtersInvisible; // Make daughters invsibile.
  G4Colour    fColour;
  LineStyle   fLineStyle;
  G4double    fLineWidth;          // Units of "normal" device linewidth, e.g.,
                                   // pixels for screen, 0.1 mm for paper.
  G4bool      fForceDrawingStyle;  // To override view parameters.
  ForcedDrawingStyle fForcedStyle;
  // Available attributes - end snippet.
};

#include "G4VisAttributes.icc"

#endif
