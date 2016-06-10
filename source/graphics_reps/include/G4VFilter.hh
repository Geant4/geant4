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
// $Id: G4VFilter.hh 66376 2012-12-18 09:42:59Z gcosmo $
//
// Abstract filter class.
//
// Jane Tinslay, March 2006
//
#ifndef G4VFILTER_HH
#define G4VFILTER_HH

#include "globals.hh"
#include "G4String.hh"
#include <iostream>

template <typename T>
class G4VFilter {

public: // With description

  typedef T Type;

  // Construct with filter name
  G4VFilter(const G4String& name);

  virtual ~G4VFilter();
  
  // Filter method
  virtual G4bool Accept(const T&) const = 0;

  // Print configuration
  virtual void PrintAll(std::ostream& ostr) const = 0;
  
  // Reset 
  virtual void Reset() = 0;

  // Filter name
  G4String Name() const;
  G4String GetName() const;

private:

  // Data member
  G4String fName;

};

template <typename T>
G4VFilter<T>::G4VFilter(const G4String& name)
  :fName(name) 
{}

template <typename T>
G4VFilter<T>::~G4VFilter() {}

template <typename T>
G4String 
G4VFilter<T>::Name() const 
{
  return fName;
}

template <typename T>
G4String 
G4VFilter<T>::GetName() const 
{
  return Name();
}

#endif

