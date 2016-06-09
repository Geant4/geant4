//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Colour.cc,v 1.6 2003/09/18 11:03:12 johna Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
// John Allison 20th October 1996

#include "G4Colour.hh"
#include "G4ios.hh"

G4Colour::G4Colour (G4double r, G4double g, G4double b, G4double a):
red (r), green (g), blue (b), alpha (a)
{
  if( red   > 1.0 ){red   = 1.0;} if( red   < 0.0 ){red   = 0.0;}
  if( green > 1.0 ){green = 1.0;} if( green < 0.0 ){green = 0.0;}
  if( blue  > 1.0 ){blue  = 1.0;} if( blue  < 0.0 ){blue  = 0.0;}
  if( alpha > 1.0 ){alpha = 1.0;} if( alpha < 0.0 ){alpha = 0.0;}
}

G4Colour::G4Colour (G4ThreeVector v):
red (v.x()), green (v.y()), blue (v.z()), alpha (1.)
{
  if( red   > 1.0 ){red   = 1.0;} if( red   < 0.0 ){red   = 0.0;}
  if( green > 1.0 ){green = 1.0;} if( green < 0.0 ){green = 0.0;}
  if( blue  > 1.0 ){blue  = 1.0;} if( blue  < 0.0 ){blue  = 0.0;}
}

G4Colour::operator G4ThreeVector() {
  return G4ThreeVector(red,green,blue);
}

std::ostream& operator << (std::ostream& os, const G4Colour& c) {
  return os << '(' << c.red << ',' << c.green << ',' << c.blue
	    << ',' << c.alpha << ')';
}

G4bool G4Colour::operator != (const G4Colour& c) const {
  if (
      (red   != c.red)   ||
      (green != c.green) ||
      (blue  != c.blue)  ||
      (alpha != c.alpha)
      )
    return true;
  return false;
}
