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
// $Id: G4VModelFactory.hh,v 1.1 2005/11/21 05:44:44 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

#endif

