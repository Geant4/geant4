// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Colour.hh,v 1.2 1999-05-25 09:10:07 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison 20th October 1996

#ifndef G4COLOUR_HH
#define G4COLOUR_HH

#include "globals.hh"
class ostream;

class G4Colour {
  friend ostream& operator << (ostream& os, const G4Colour& c);
public:
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
