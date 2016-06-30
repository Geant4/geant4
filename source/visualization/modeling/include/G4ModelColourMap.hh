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
// $Id: G4ModelColourMap.hh 95593 2016-02-16 10:48:50Z gcosmo $
//
// Generic variable->G4Colour map, where "variable" is the template
// parameter.
//
// Jane Tinslay March 2006

#ifndef G4MODELCOLOURMAP_HH
#define G4MODELCOLOURMAP_HH

#include "G4Colour.hh"
#include "G4String.hh"
#include <map>

template <typename T>
class G4ModelColourMap {

public: // With description

  G4ModelColourMap();

  virtual ~G4ModelColourMap();

  // Configuration functions
  void Set(const T&, const G4Colour&);
  void Set(const T&, const G4String&);
  G4Colour& operator[](const T& quantity);

  // Access functions
  const std::map<T, G4Colour>& GetBasicMap() const;
  bool GetColour(const T&, G4Colour&) const;
  void Print(std::ostream& ostr) const;

private:

  // Data member
  std::map<T, G4Colour> fMap;

};

template <typename T>
G4Colour& 
G4ModelColourMap<T>::operator[](const T& quantity) {return fMap[quantity];}

template <typename T>
G4ModelColourMap<T>::G4ModelColourMap() {}

template <typename T>
G4ModelColourMap<T>::~G4ModelColourMap() {}

template <typename T>
void
G4ModelColourMap<T>::Set(const T& quantity, const G4String& colour) 
{
  G4Colour myColour;
  
  // Will not setup the map if colour key does not exist
  if (!G4Colour::GetColour(colour, myColour)) {
    G4ExceptionDescription ed;
    ed << "G4Colour with key "<<colour<<" does not exist ";
    G4Exception
      ("G4ColourMap::Set(Charge charge, const G4String& colour)",
       "modeling0108", JustWarning, ed);
    return;
  }
  
  
  // Will not modify myColour if colour key does not exist
  Set(quantity, myColour);
}

template <typename T>
void
G4ModelColourMap<T>::Set(const T& quantity, const G4Colour& colour) 
{
  fMap[quantity] = colour;
}   

template <typename T>
const std::map<T, G4Colour>& G4ModelColourMap<T>::GetBasicMap() const
{ return fMap; }

template <typename T>
bool
G4ModelColourMap<T>::GetColour(const T& quantity, G4Colour& colour) const
{
  typename std::map<T, G4Colour>::const_iterator iter = fMap.find(quantity);

  if (iter != fMap.end()) {
    colour = iter->second;
    return true;
  }
  
  return false;
}

template <typename T>
void
G4ModelColourMap<T>::Print(std::ostream& ostr) const
{
  typename std::map<T, G4Colour>::const_iterator iter = fMap.begin();
  
  while (iter != fMap.end()) {
    ostr<< iter->first <<" : "<< iter->second <<G4endl;
    iter++;
  }
}
#endif
