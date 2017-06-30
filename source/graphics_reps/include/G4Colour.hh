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
// $Id: G4Colour.hh 104093 2017-05-11 07:29:45Z gcosmo $
//
// 
// John Allison 20th October 1996

// Class Description:
// Class G4Colour has 4 fields, which represent the RGBA (red, green, blue, 
// and alpha) components of colour. Each component takes a value between 
// 0 and 1. If an irrelevant value, i.e., a value less than 0 or greater 
// than 1, is given as an argument of the constructor, such a value is 
// automatically clipped to 0 or 1. Alpha is opacity (1 = opaque).
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
//      G4Colour  brown   (0.45,0.25,0.0) ;  // G4 logo brown
//      G4Colour  red     (1.0, 0.0, 0.0) ;  // red
//      G4Colour  green   (0.0, 1.0, 0.0) ;  // green
//      G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
//      G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
//      G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
//      G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
// 
// For convenience, static member functions are also defined for the above colours. 
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
#include "G4ThreeVector.hh"
#include <iostream>
#include <map>

class G4Colour {

  friend std::ostream& operator << (std::ostream&, const G4Colour&);

public: // With description

  G4Colour (G4double r = 1., G4double g = 1., G4double b = 1.,
            G4double a = 1.);

  G4Colour (G4ThreeVector);
  // Converts the components of the 3-vector into red, green, blue.
  // The opacity, alpha = 1.

  operator G4ThreeVector();
  // Converts red, green, blue into the components of a 3-vector.

  G4bool operator != (const G4Colour& c) const;
  G4bool operator == (const G4Colour& c) const {return !(operator != (c));}

  G4Colour& operator += (const G4Colour& rhs) {*this = rhs; return *this;}
  // Note: This is required by RayTracer in its use of G4THitsMap.
  // Adding colours, without also taking brightness into account, does not make
  // sense, so let us make it synonymous with operator=, which is, I guess,
  // equivalent to covering the old colour with the new, like a coat of paint.

  G4double GetRed   () const;
  G4double GetGreen () const;
  G4double GetBlue  () const;
  G4double GetAlpha () const;  // alpha = opacity = 1. - transparency.

  void SetRed   (G4double);
  void SetGreen (G4double);
  void SetBlue  (G4double);
  void SetAlpha (G4double);  // alpha = opacity = 1. - transparency.

  static G4Colour White();
  static G4Colour Gray();
  static G4Colour Grey();
  static G4Colour Black();
  static G4Colour Brown();  // G4 logo brown
  static G4Colour Red();
  static G4Colour Green();
  static G4Colour Blue(); 
  static G4Colour Cyan();
  static G4Colour Magenta();
  static G4Colour Yellow();

  static void AddToMap(const G4String& key, const G4Colour& colour);
  // Add user defined colour to colour map with given key. Standard
  // colours are added to map by default.

  static G4bool GetColour(const G4String& key, G4Colour& result);
  // Get colour for given key. Returns false if key doesn't exist 
  // in colour map, leaving result unchanged. Colour map
  // is not sensitive to key case.

  static const std::map<G4String, G4Colour>& GetMap();
  
private:
  G4double red, green, blue, alpha;

  static G4ThreadLocal std::map<G4String, G4Colour> *fColourMap;
  static G4ThreadLocal G4bool fInitColourMap;
  static void InitialiseColourMap();
    
};

inline G4double G4Colour::GetRed   () const {return red;}
inline G4double G4Colour::GetGreen () const {return green;}
inline G4double G4Colour::GetBlue  () const {return blue;}
inline G4double G4Colour::GetAlpha () const {return alpha;}
inline G4Colour G4Colour::White()   {return G4Colour(1.0, 1.0, 1.0);}
inline G4Colour G4Colour::Gray()    {return G4Colour(0.5, 0.5, 0.5);}
inline G4Colour G4Colour::Grey()    {return G4Colour(0.5, 0.5, 0.5);}
inline G4Colour G4Colour::Black()   {return G4Colour(0.0, 0.0, 0.0);}
inline G4Colour G4Colour::Brown()   {return G4Colour(0.45,0.25,0.0);}
inline G4Colour G4Colour::Red()     {return G4Colour(1.0, 0.0, 0.0);}
inline G4Colour G4Colour::Green()   {return G4Colour(0.0, 1.0, 0.0);}
inline G4Colour G4Colour::Blue()    {return G4Colour(0.0, 0.0, 1.0);} 
inline G4Colour G4Colour::Cyan()    {return G4Colour(0.0, 1.0, 1.0);}
inline G4Colour G4Colour::Magenta() {return G4Colour(1.0, 0.0, 1.0);}
inline G4Colour G4Colour::Yellow()  {return G4Colour(1.0, 1.0, 0.0);}

#endif
