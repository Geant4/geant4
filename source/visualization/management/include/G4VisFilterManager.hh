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
// $Id: G4VisFilterManager.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Filter manager. Manages filter models, factories, messengers, 
// command placement, filter mode etc
//
// Jane Tinslay, March 2006
//
#ifndef G4VISFILTERMANAGER_HH
#define G4VISFILTERMANAGER_HH

#include "G4String.hh"
#include "G4UImessenger.hh"
#include "G4VFilter.hh"
#include "G4VModelFactory.hh"
#include <vector>

namespace FilterMode {
  enum Mode {Soft, Hard};
}

template <typename T>
class G4VisFilterManager {

public:

  // Construct with command placement
  G4VisFilterManager(const G4String&);

  virtual ~G4VisFilterManager();

  // Useful typedef's
  typedef G4VFilter<T> Filter;
  typedef G4VModelFactory<Filter> Factory;

  // Registration methods
  void Register(Filter*);
  void Register(Factory*); 

  // Do filtering
  bool Accept(const T&);

  // Command placement
  G4String Placement() const;

  // Filter mode operations
  void SetMode(const FilterMode::Mode&);
  void SetMode(const G4String&);
  FilterMode::Mode GetMode() const;

  // Print configuration
  void Print(std::ostream& ostr, const G4String& name="") const;

  // Accessors
  const std::vector<Filter*>& FilterList() const;
  const std::vector<Factory*>& FactoryList() const;

private:

  // Data members
  G4String fPlacement; // Placement 
  FilterMode::Mode fMode;
  std::vector<Factory*> fFactoryList;
  std::vector<Filter*> fFilterList;
  std::vector<G4UImessenger*> fMessengerList;

};

template <typename T>
G4VisFilterManager<T>::G4VisFilterManager(const G4String& placement)
  :fPlacement(placement)
{
  fMode = FilterMode::Hard;
}

template <typename T>
G4VisFilterManager<T>::~G4VisFilterManager() 
{
  // Cleanup
  std::vector<G4UImessenger*>::iterator iterMsgr = fMessengerList.begin();
  
  while (iterMsgr != fMessengerList.end()) {
    delete *iterMsgr;
    iterMsgr++;
  }
  
  typename std::vector<Factory*>::iterator iterFactory = fFactoryList.begin();
  
  while (iterFactory != fFactoryList.end()) {
    delete *iterFactory;       
    iterFactory++;
  }

  typename std::vector<Filter*>::iterator iterFilter = fFilterList.begin();
  
  while (iterFilter != fFilterList.end()) {
    delete *iterFilter;       
    iterFilter++;
  }
}

template <typename T>
void
G4VisFilterManager<T>::Register(Filter* filter)
{
  fFilterList.push_back(filter);
}

template <typename T>
void
G4VisFilterManager<T>::Register(Factory* factory)
{
  fFactoryList.push_back(factory);

  fMessengerList.push_back(new G4VisCommandModelCreate<Factory>(factory, fPlacement));
}

template <typename T>
bool
G4VisFilterManager<T>::Accept(const T& obj)
{
  typename std::vector<Filter*>::const_iterator iter = fFilterList.begin();
  bool passed(true);
  
  while (passed && (iter != fFilterList.end())) {
    passed = (*iter)->Accept(obj);
    iter++;
  }

  return passed;
}

template <typename T>
G4String
G4VisFilterManager<T>::Placement() const
{
  return fPlacement;
}

template <typename T>
void
G4VisFilterManager<T>::SetMode(const G4String& mode) 
{
  bool result(false);
  
  G4String myMode(mode);
  myMode.toLower();

  if (myMode == "soft") {result = true; SetMode(FilterMode::Soft);}
  else if (myMode == "hard") {result = true; SetMode(FilterMode::Hard);}

  if (!result) {
    G4ExceptionDescription ed;
    ed << "Invalid Filter mode: "<<mode;
    G4Exception
      ("G4VisFilterManager::SetMode(const G4String& mode)", "visman0101", JustWarning, ed);
  }
}

template <typename T>
void
G4VisFilterManager<T>::SetMode(const FilterMode::Mode& mode) 
{
  fMode = mode;
}

template <typename T>
FilterMode::Mode
G4VisFilterManager<T>::GetMode() const
{
  return fMode;
}

template <typename T>
void
G4VisFilterManager<T>::Print(std::ostream& ostr, const G4String& name) const
{ 
  ostr<<"Registered filter factories:"<<std::endl;
  typename std::vector<Factory*>::const_iterator iterFactory = fFactoryList.begin();

  while (iterFactory != fFactoryList.end()) {
    (*iterFactory)->Print(ostr);
    iterFactory++;
  }

  if (0 == fFactoryList.size()) ostr<<"  None"<<std::endl;

  ostr<<std::endl;
  ostr<<"Registered filters:"<<std::endl;

  typename std::vector<Filter*>::const_iterator iterFilter = fFilterList.begin();

  while (iterFilter != fFilterList.end()) {
    if (!name.isNull()) {
      if ((*iterFilter)->Name() == name) (*iterFilter)->PrintAll(ostr);
    }
    else {
      (*iterFilter)->PrintAll(ostr);
    }
    iterFilter++;
  }

  if (0 == fFilterList.size()) ostr<<"  None"<<std::endl;
}

template <typename T>
const std::vector< G4VFilter<T>* >&
G4VisFilterManager<T>::FilterList() const
{ 
  return fFilterList;
}

template <typename T>
const std::vector< G4VModelFactory< G4VFilter<T> >* >&
G4VisFilterManager<T>::FactoryList() const
{ 
  return fFactoryList;
}

#endif
