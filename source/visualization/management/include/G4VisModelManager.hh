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
// $Id: G4VisModelManager.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Generic model manager. Manages models, associated
// factories, messengers, command placement etc
//
// Jane Tinslay, March 2006
//
#ifndef G4VISMODELMANAGER_HH
#define G4VISMODELMANAGER_HH

#include "G4UImessenger.hh"
#include "G4VisCommandModelCreate.hh"
#include "G4VisListManager.hh"
#include "G4VModelFactory.hh"
#include <vector>

template <typename Model>
class G4VisModelManager {

public: // With description

  // Useful typedef's
  typedef G4VisListManager<Model> List;
  typedef G4VModelFactory<Model> Factory;

  G4VisModelManager(const G4String&);
  virtual ~G4VisModelManager();

  // Registration methods
  void Register(Model*);
  void Register(Factory*); 

  // Change/Retrieve "Current" object
  void SetCurrent(const G4String&);
  const Model* Current() const;

  // Command placement
  G4String Placement() const;

  // Print factory and model data
  void Print(std::ostream& ostr, const G4String& name="") const;

  // Accessors
  const List* ListManager() const;
  const std::vector<Factory*>& FactoryList() const;

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps Coverity happy.
  G4VisModelManager (const G4VisModelManager&);
  G4VisModelManager& operator = (const G4VisModelManager&);

  // Data members
  G4String fPlacement;
  List* fpModelList;  
  std::vector<Factory*> fFactoryList;
  std::vector<G4UImessenger*> fMessengerList;

};

template <typename Model>
G4VisModelManager<Model>::G4VisModelManager(const G4String& placement)
  :fPlacement(placement)
  ,fpModelList(new List)
{}

template <typename Model>
G4VisModelManager<Model>::~G4VisModelManager() 
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

  delete fpModelList;
}

template <typename Model>
void
G4VisModelManager<Model>::Register(Model* model)
{
  fpModelList->Register(model);
}

template <typename Model>
void
G4VisModelManager<Model>::Register(Factory* factory)
{
  // Assume ownership
  fFactoryList.push_back(factory);

  // Generate "create" command for this factory
  fMessengerList.push_back(new G4VisCommandModelCreate<Factory>(factory, fPlacement));
}

template <typename Model>
void
G4VisModelManager<Model>::SetCurrent(const G4String& model) 
{
  fpModelList->SetCurrent(model);
}

template <typename Model>
const Model*
G4VisModelManager<Model>::Current() const
{
  return fpModelList->Current();
}

template <typename Model>
G4String
G4VisModelManager<Model>::Placement() const
{
  return fPlacement;
}

template <typename Model>
void
G4VisModelManager<Model>::Print(std::ostream& ostr, const G4String& name) const
{
  ostr<<"Registered model factories:"<<std::endl;

  typename std::vector<Factory*>::const_iterator iter = fFactoryList.begin();

  while (iter != fFactoryList.end()) {
    (*iter)->Print(ostr);
    iter++;
  }

  if (0 == fFactoryList.size()) ostr<<"  None"<<std::endl;

  ostr<<std::endl;
  ostr<<"Registered models: "<<std::endl;

  fpModelList->Print(ostr, name);
}

template <typename Model>
const G4VisListManager<Model>*
G4VisModelManager<Model>::ListManager() const
{
  return fpModelList;
}

template <typename Model>
const std::vector<G4VModelFactory<Model>*>&
G4VisModelManager<Model>::FactoryList() const
{
  return fFactoryList;
}

#endif
