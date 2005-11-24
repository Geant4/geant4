//
// File name:     RadmonPhysicsList.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsList.cc,v 1.3 2005-11-24 02:38:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonPhysicsList.hh"

#include "RadmonVPhysicsLayout.hh"
#include "RadmonVSubPhysicsList.hh"
#include "RadmonVSubPhysicsListFactory.hh"
#include "RadmonPhysicsInfoList.hh"

                                                RadmonPhysicsList :: RadmonPhysicsList(RadmonVPhysicsLayout * layout, RadmonVSubPhysicsListFactory * factory)
:
 physicsLayout(layout),
 subPhysicsListFactory(factory),
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
 
 AddTransportation();

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
  (*i)->ConstructProcess();
  
  i++;
 }
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
