//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4Colour.cc,v 1.11 2007-01-05 14:06:28 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison 20th October 1996

#include "G4Colour.hh"

#include <sstream>

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

map<G4String, G4Colour> G4Colour::fColourMap;
bool G4Colour::fInitColourMap = false;

void
G4Colour::AddToMap(const G4String& key, const G4Colour& colour) 
{
  // Convert to lower case since colour map is case insensitive
  G4String myKey(key);
  myKey.toLower();

  map<G4String, G4Colour>::iterator iter = fColourMap.find(myKey);
  
  if (iter == fColourMap.end()) fColourMap[myKey] = colour;  
  else {
    std::ostringstream o; 
    o << "G4Colour with key "<<myKey<<" already exists ";
    G4Exception
      ("G4Colour::AddToMap(const G4String& key, const G4Colour& colour)",
       "ColourKeyExists", JustWarning, o.str().c_str());
  }
}

void
G4Colour::InitialiseColourMap() 
{
  // Standard colours
  AddToMap("white",   G4Colour::White());
  AddToMap("gray",    G4Colour::Gray());
  AddToMap("grey",    G4Colour::Grey());
  AddToMap("black",   G4Colour::Black());
  AddToMap("red",     G4Colour::Red());
  AddToMap("green",   G4Colour::Green());
  AddToMap("blue",    G4Colour::Blue());
  AddToMap("cyan",    G4Colour::Cyan());
  AddToMap("magenta", G4Colour::Magenta());
  AddToMap("yellow",  G4Colour::Yellow());
}

bool
G4Colour::GetColour(const G4String& key, G4Colour& result) 
{
  if (false == fInitColourMap) {
    fInitColourMap = true;
    // Add standard colours to map
    InitialiseColourMap();
  }
 
  G4String myKey(key);
  myKey.toLower();
 
  map<G4String, G4Colour>::iterator iter = fColourMap.find(myKey);

  // Don't modify "result" if colour was not found in map
  if (iter == fColourMap.end()) return false;
  
  result = iter->second;

  return true;
}
