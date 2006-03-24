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
// $Id: G4VisModelManager.hh,v 1.1 2006-03-24 22:03:02 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  void Print(std::ostream& ostr, const G4String& name) const;

private:

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
  fpModelList->Print(ostr, name);
}

#endif
