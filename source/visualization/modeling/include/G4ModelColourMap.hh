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
// $Id: G4ModelColourMap.hh,v 1.1 2006-03-17 03:24:02 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    std::ostringstream o;
    o << "G4Colour with key "<<colour<<" does not exist ";
    G4Exception
      ("G4ColourMap::Set(Charge charge, const G4String& colour)",
       "NonExistentColour", JustWarning, o.str().c_str());
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
