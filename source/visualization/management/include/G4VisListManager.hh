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
// $Id: G4VisListManager.hh,v 1.3 2005/11/23 20:25:22 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
#include <sstream>

template <typename T>
class G4VisListManager {

public: // With description
  
  G4VisListManager();

  virtual ~G4VisListManager();

  void Register(T* ptr);
  // Register ptr. Manager assumes ownership and
  // ptr becomes current

  void SetCurrent(const G4String& name);
  const T* Current() const {return fpCurrent;}
  
  void Print(std::ostream& ostr, const G4String& name) const;

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

  typename std::map<G4String, T*>::const_iterator iter = fMap.find(ptr->Name());
 
  if (iter == fMap.end()) {
    fMap[ptr->Name()] = ptr;
    fpCurrent = ptr;    
  }
  else {
    std::ostringstream o;
    o << "Key "<<ptr->Name()<<" already registered";
    G4Exception
      ("G4VisListManager<T>::Register(T* ptr) ",
       "KeyExists", FatalErrorInArgument, o.str().c_str());
  }
}

template <typename T>
void
G4VisListManager<T>::SetCurrent(const G4String& name)
{
  typename std::map<G4String, T*>::const_iterator iter = fMap.find(name);

  if (iter != fMap.end()) fpCurrent = fMap[name];
  else {
    std::ostringstream o;
    o << "Key "<<name<<" has not been registered";
    G4Exception
      ("G4VisListManager<T>::SetCurrent(T* ptr) ",
       "NonExistentName", FatalErrorInArgument, o.str().c_str());
  }
}

template <typename T>
void
G4VisListManager<T>::Print(std::ostream& ostr, const G4String& name) const
{
  ostr<<"Current: "<<fpCurrent->Name()<<std::endl;

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
      ostr<<G4endl;
      iter++;
    }
  }
}

#endif
