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
//
// File name:     RadmonPhysicsList.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsList.cc,v 1.6 2006/06/29 16:18:57 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//

// Include files
#include "RadmonPhysicsList.hh"

#include "RadmonVPhysicsLayout.hh"
#include "RadmonVSubPhysicsList.hh"
#include "RadmonVSubPhysicsListFactory.hh"
#include "RadmonPhysicsInfoList.hh"
#include "RadmonSteppingAction.hh"
#include "RadmonPhysicsSteppingAction.hh"

                                                RadmonPhysicsList :: RadmonPhysicsList(RadmonVPhysicsLayout * layout, RadmonVSubPhysicsListFactory * factory)
:
 physicsLayout(layout),
 subPhysicsListFactory(factory),
 steppingAction(0),
 initializationMethodsCalled(false),
 changed(true)
{
 if (physicsLayout==0)
  G4Exception("RadmonPhysicsList::RadmonPhysicsList: layout==0.");

 if (factory==0)
  G4Exception("RadmonPhysicsList::RadmonPhysicsList: factory==0.");
 
 physicsLayout->AttachObserver(this);
}



                                                RadmonPhysicsList :: ~RadmonPhysicsList()
{
 physicsLayout->DetachObserver(this);
 
 if (steppingAction)
 {
  RadmonSteppingAction::Instance()->DetachObserver(steppingAction);
  delete steppingAction;
 }
 
 Destruct();
 
 delete subPhysicsListFactory;
}





void                                            RadmonPhysicsList :: OnLayoutChange(void)
{
 if (initializationMethodsCalled)
 {
  G4cout << "RadmonPhysicsList::OnLayoutChange: Physics initialization done, changes to the physics will be ignored" << G4endl;
  return;
 }
 
 changed=true;
}





void                                            RadmonPhysicsList :: ConstructParticle(void)
{
}



void                                            RadmonPhysicsList :: ConstructProcess(void)
{
 CheckUpdate();
 initializationMethodsCalled=true;
 
 if (!steppingAction)
 {
  steppingAction=new RadmonPhysicsSteppingAction;
  RadmonSteppingAction::Instance()->AttachObserver(steppingAction);
 }
 
 SubPhysiscsLists::iterator i(subPhysiscsLists.begin());
 const SubPhysiscsLists::iterator end(subPhysiscsLists.end());
 
 while (i!=end)
 {
  (*i)->ConstructParticle();
  
  i++;
 }
 
 UpdateProcessManagers();

 i=subPhysiscsLists.begin();
 while (i!=end)
 {
  (*i)->ConstructProcess();
  UpdateProcessManagers();
  
  i++;
 }

 AddTransportation();
}



void                                            RadmonPhysicsList :: SetCuts(void)
{
 CheckUpdate();
 initializationMethodsCalled=true;

 SubPhysiscsLists::iterator i(subPhysiscsLists.begin());
 const SubPhysiscsLists::iterator end(subPhysiscsLists.end());
 
 while (i!=end)
 {
  (*i)->ConstructParticle();
  
  i++;
 }

 i=subPhysiscsLists.begin();
 while (i!=end)
 {
  (*i)->SetCuts();
  
  i++;
 }
}





void                                            RadmonPhysicsList :: Destruct(void)
{
 SubPhysiscsLists::iterator i(subPhysiscsLists.begin());
 const SubPhysiscsLists::iterator end(subPhysiscsLists.end());
 
 while (i!=end)
 {
  delete (*i);
  
  i++;
 }
 
 subPhysiscsLists.clear();
}





void                                            RadmonPhysicsList :: CheckUpdate(void)
{
 if (!changed)
  return;
  
 changed=false;
 
 Destruct();
 
 const G4int n(physicsLayout->GetNPhysicsLists());
 
 for (G4int i(0); i<n; i++)
 {
  G4String name(physicsLayout->GetPhysicsListName(i));
  
  RadmonVSubPhysicsList * subPhysicsList(subPhysicsListFactory->CreateSubPhysicsList(name));
  
  if (subPhysicsList!=0)
  {
   G4int j(physicsLayout->GetPhysicsListNAttributes(name));
   
   while (j>0)
   {
    j--;
    
    const G4String & attributeName(physicsLayout->GetPhysicsListAttributeName(name, j));
    subPhysicsList->SetPhysicsListAttribute(attributeName, physicsLayout->GetPhysicsListAttribute(name, attributeName));
   }
   
   SubPhysiscsLists::const_iterator k(subPhysiscsLists.begin());
   const SubPhysiscsLists::const_iterator end(subPhysiscsLists.end());
   
   while (k!=end)
   {
    if (subPhysicsList->Provides().CollidesWith((*k)->Provides()))
     G4cout << "RadmonPhysicsList::OnLayoutChange: Warning physics list \"" << name << "\" provides features common to other physics list." << G4endl;

    k++;
   }
   
   subPhysiscsLists.push_back(subPhysicsList);
  }
  else
   G4cout << "RadmonPhysicsList::OnLayoutChange: Physics list with name \"" << name << "\" not found." << G4endl;
 }
}





void                                            RadmonPhysicsList :: UpdateProcessManagers(void)
{
 G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
 G4ParticleTable::G4PTblDicIterator & particleIterator(*particleTable->GetIterator());
 
 particleIterator.reset();
 while (particleIterator())
 {
  G4ParticleDefinition * particle(particleIterator.value());
  
  if (particle->GetProcessManager()==0)
   AddProcessManager(particle);
 }
}

