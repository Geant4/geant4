// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Colour.hh,v 1.4 1999-12-15 14:50:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison 20th October 1996

// Class Description:
// Class G4Colour has 4 fields, which represent the RGBA (red, green, blue, 
// and alpha) components of colour. Each component takes a value between 
// 0 and 1. If an irrelevant value, i.e., a value less than 0 or greater 
// than 1, is given as an argument of the constructor, such a value is 
// automatically clipped to 0 or 1. Alpha is opacity, which is not used 
// at present. You can use its default value 1, which means "opaque" 
// in instantiation of G4Colour. 
// 
// A G4Colour object is instantiated by giving red, green, and blue 
// components to its constructor, i.e., 
// 
//      G4Colour::G4Colour ( G4double r = 1.0, 
//                           G4double g = 1.0, 
//                           G4double b = 1.0, 
//                           G4double a = 1.0);
//                                   // 0<=red, green, blue <= 1.0 
// 
// The default value of each component is 1.0. That is to say, the default 
// colour is "white".  For example, colours which are often used can be 
// instantiated as follows: 
// 
//      G4Colour  white   ()              ;  // white
//      G4Colour  white   (1.0, 1.0, 1.0) ;  // white
//      G4Colour  gray    (0.5, 0.5, 0.5) ;  // gray
//      G4Colour  black   (0.0, 0.0, 0.0) ;  // black
//      G4Colour  red     (1.0, 0.0, 0.0) ;  // red
//      G4Colour  green   (0.0, 1.0, 0.0) ;  // green
//      G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
//      G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
//      G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
//      G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
// 
// After instantiation of a G4Colour object, you can access to its components 
// with the following access functions: 
// 
//      G4double G4Colour::GetRed   () const ; // Get the red   component. 
//      G4double G4Colour::GetGreen () const ; // Get the green component.
//      G4double G4Colour::GetBlue  () const ; // Get the blue  component.
//
// Class Description - End:

#ifndef G4COLOUR_HH
#define G4COLOUR_HH

#include "globals.hh"
class G4std::ostream;

class G4Colour {
  friend G4std::ostream& operator << (G4std::ostream& os, const G4Colour& c);

public: // With description

  G4Colour (G4double r = 1., G4double g = 1., G4double b = 1.,
	    G4double a = 1.);
  G4bool operator != (const G4Colour& c) const;
  G4double GetRed   () const;
  G4double GetGreen () const;
  G4double GetBlue  () const;
  G4double GetAlpha () const;  // alpha = opacity = 1. - transparency.
private:
  G4double red, green, blue, alpha;
};

inline G4Colour::G4Colour (G4double r, G4double g, G4double b, G4double a):
red (r), green (g), blue (b), alpha (a)
{
  if( red   > 1.0 ){red   = 1.0;} if( red   < 0.0 ){red   = 0.0;}
  if( green > 1.0 ){green = 1.0;} if( green < 0.0 ){green = 0.0;}
  if( blue  > 1.0 ){blue  = 1.0;} if( blue  < 0.0 ){blue  = 0.0;}
  if( alpha > 1.0 ){alpha = 1.0;} if( alpha < 0.0 ){alpha = 0.0;}
}

inline G4double G4Colour::GetRed   () const {return red;}
inline G4double G4Colour::GetGreen () const {return green;}
inline G4double G4Colour::GetBlue  () const {return blue;}
inline G4double G4Colour::GetAlpha () const {return alpha;}

#endif
