// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisAttributes.hh,v 1.6 2000-05-13 09:33:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  23rd October 1996

// Class Description:
// Visualization attributes are a set of information associated with the 
// visualizable objects. This information is necessary only for 
// visualization, and is not included in geometrical information such 
// as shapes, position, and orientation. 
// A typical example of a visualization attribute is "colour". 
// For example, in visualizing a box, the Visualization Manager must know 
// its colour. If an object to be visualized has not been assigned a set of 
// visualization attributes, then a proper default set is used 
// automatically. A set of visualization attributes is held by an 
// instance of class G4VisAttributes defined in the graphics_reps 
// category. The followings are commonly-used attributes:
//   - visibility
//   - force wireframe style, force solid style
//   - colour 
// Class Description - End:


#ifndef __G4VISATTRIBUTES_HH__
#define __G4VISATTRIBUTES_HH__

#include "G4ios.hh"
#include "globals.hh"
#include "G4Colour.hh"
#include "G4Color.hh"

class G4VisAttributes {

  friend G4std::ostream& operator << (G4std::ostream& os, const G4VisAttributes& a);

public: // With description

  enum LineStyle {unbroken, dashed, dotted};
  enum ForcedDrawingStyle {wireframe, solid};

  G4VisAttributes ();
  G4VisAttributes (G4bool visibility);
  G4VisAttributes (const G4Colour& colour);
  G4VisAttributes (G4bool visibility, const G4Colour& colour);

  static const G4VisAttributes Invisible;

  G4bool operator != (const G4VisAttributes& a) const;
  G4bool operator == (const G4VisAttributes& a) const;

  G4bool          IsVisible                      () const;
  G4bool          IsDaughtersInvisible           () const;
  const G4Colour& GetColour                      () const;
  const G4Color&  GetColor                       () const;
  LineStyle       GetLineStyle                   () const;
  G4double        GetLineWidth                   () const;
  G4bool          IsForceDrawingStyle            () const;
  ForcedDrawingStyle GetForcedDrawingStyle () const;

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

private:

  G4bool      fVisible;            // Visibility flag
  G4bool      fDaughtersInvisible; // Make daughters invsibile.
  G4Colour    fColour;
  LineStyle   fLineStyle;
  G4double    fLineWidth;          // Units of "normal" device linewidth, e.g.,
                                   // pixels for screen, 0.1 mm for paper.
  G4bool      fForceDrawingStyle;  // To override view parameters.
  ForcedDrawingStyle fForcedStyle;
};

#include "G4VisAttributes.icc"

#endif
