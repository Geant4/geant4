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
// $Id: G4VModelFactory.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Abstract base class for model factories. Derived classes
// must implement the Create method, creating a new model
// and associated messengers.
// Class Description - End:

#ifndef G4VMODELFACTORY
#define G4VMODELFACTORY

#include "G4String.hh"
#include <iostream>
#include <utility>
#include <vector>

class G4UImessenger;

template <class T>
class G4VModelFactory {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair<T*, Messengers> ModelAndMessengers;
  
  G4VModelFactory(const G4String& name);

  virtual ~G4VModelFactory();

  G4String Name();

  virtual ModelAndMessengers Create(const G4String& placement, const G4String& modelName) = 0;

  void Print(std::ostream& ostr) const;

private:

  G4String fName;

};

template <typename T>
G4VModelFactory<T>::G4VModelFactory(const G4String& name)
  :fName(name) 
{}

template <typename T>
G4VModelFactory<T>::~G4VModelFactory() {}

template <typename T>
G4String 
G4VModelFactory<T>::Name() 
{
  return fName;
}

template <typename T>
void 
G4VModelFactory<T>::Print(std::ostream& ostr) const
{
  ostr<<"  "<<fName<<std::endl; 
} 

#endif

