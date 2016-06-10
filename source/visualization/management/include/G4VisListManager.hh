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
// $Id: G4VisListManager.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Templated class to manage a set of objects, with one named
// as current. Owns registered objects.
// Class Description - End:

#ifndef G4VISLISTMANAGER_HH
#define G4VISLISTMANAGER_HH

#include "G4String.hh"
#include <map>
#include <ostream>

template <typename T>
class G4VisListManager {

public: // With description
  
  G4VisListManager();

  virtual ~G4VisListManager();

  // Register ptr. Manager assumes ownership and
  // ptr becomes current
  void Register(T* ptr);

  void SetCurrent(const G4String& name);

  // Accessors
  const T* Current() const {return fpCurrent;}
  const std::map<G4String, T*>& Map() const;

  // Print configuration
  void Print(std::ostream& ostr, const G4String& name="") const;

private:

  // Data members
  std::map<G4String, T*> fMap;
  T* fpCurrent;

};

template <typename T>
G4VisListManager<T>::G4VisListManager()
  :fpCurrent(0) 
{}

template <typename T>
G4VisListManager<T>::~G4VisListManager() 
{
  typename std::map<G4String, T*>::iterator iter = fMap.begin();

  while (iter != fMap.end()) {
    delete iter->second;    
    iter++;
  }
}

template <typename T>
void
G4VisListManager<T>::Register(T* ptr) 
{
  assert (0 != ptr);

  // Add to map.  Replace if name the same.
  fMap[ptr->Name()] = ptr;
  fpCurrent = ptr;    
}

template <typename T>
void
G4VisListManager<T>::SetCurrent(const G4String& name)
{
  typename std::map<G4String, T*>::const_iterator iter = fMap.find(name);

  if (iter != fMap.end()) fpCurrent = fMap[name];
  else {
    G4ExceptionDescription ed;
    ed << "Key \"" << name << "\" has not been registered";
    G4Exception
      ("G4VisListManager<T>::SetCurrent(T* ptr) ",
       "visman0102", JustWarning, ed, "Non-existent name");
  }
}

template <typename T>
void
G4VisListManager<T>::Print(std::ostream& ostr, const G4String& name) const
{
  if (0 == fMap.size()) {
    G4cout<<"  None"<<std::endl;
    return;
  }
    
  ostr<<"  Current: "<<fpCurrent->Name()<<std::endl;

  if (!name.isNull()) {
    // Print out specified object
    typename std::map<G4String, T*>::const_iterator iter = fMap.find(name);

    if (iter != fMap.end()) {
      iter->second->Print(ostr);
    }
    else {
      ostr<<name<<" not found "<<std::endl;
    }
  }
  else {
    typename std::map<G4String, T*>::const_iterator iter = fMap.begin();
    while (iter != fMap.end()) {
      iter->second->Print(ostr);
      ostr<<std::endl;
      iter++;
    }
  }
}

template <typename T>
const std::map<G4String, T*>&
G4VisListManager<T>::Map() const
{
  return fMap;
}

#endif
