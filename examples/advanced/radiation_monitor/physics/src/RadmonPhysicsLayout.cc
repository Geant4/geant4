//
// File name:     RadmonPhysicsLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsLayout.cc,v 1.1 2005-11-07 17:52:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonPhysicsLayout.hh"

                                                RadmonPhysicsLayout :: RadmonPhysicsLayout()
{
}


                                                
                                                RadmonPhysicsLayout :: ~RadmonPhysicsLayout()
{
}





void                                            RadmonPhysicsLayout :: AddPhysicsList(const G4String & physicsListName)
{
 if (subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::AddPhysicsList: Physics list \"" << physicsListName << "\" is just present." << G4endl;
  return;
 }
  
 RadmonSubPhysicsListLayout & physicsList(subPhysicsListCollection.AppendItem());
 physicsList.SetLabel(physicsListName);
 NotifyChange();
}



void                                            RadmonPhysicsLayout :: RemovePhysicsList(const G4String & physicsListName)
{
 if (!subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::RemovePhysicsList: Physics list \"" << physicsListName << "\" is not present." << G4endl;
  return;
 }
 
 subPhysicsListCollection.RemoveItemByLabel(physicsListName);
 NotifyChange();
}



G4int                                           RadmonPhysicsLayout :: GetNPhysicsLists(void) const
{
 return subPhysicsListCollection.GetNItems();
}



const G4String &                                RadmonPhysicsLayout :: GetPhysicsListName(G4int index) const
{
 return subPhysicsListCollection.GetItem(index).GetLabel();
}





void                                            RadmonPhysicsLayout :: SetPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName, const G4String & attributeValue)
{
 if (!subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::SetPhysicsListAttribute: Physics list \"" << physicsListName << "\" is not present." << G4endl;
  return;
 }
 
 RadmonSubPhysicsListLayout & subPhysicsList(subPhysicsListCollection.FindItemByLabel(physicsListName));
 
 if (subPhysicsList.GetAttribute(attributeName, attributeValue+"#")==attributeValue)
  return;

 subPhysicsList.SetAttribute(attributeName, attributeValue);
 NotifyChange();
}



void                                            RadmonPhysicsLayout :: ClearPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName)
{
 if (!subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::ClearPhysicsListAttribute: Physics list \"" << physicsListName << "\" is not present." << G4endl;
  return;
 }
 
 RadmonSubPhysicsListLayout & subPhysicsList(subPhysicsListCollection.FindItemByLabel(physicsListName));
 subPhysicsList.ClearAttribute(attributeName);
}



G4int                                           RadmonPhysicsLayout :: GetPhysicsListNAttributes(const G4String & physicsListName) const
{
 if (!subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::GetPhysicsListNAttributes: Physics list \"" << physicsListName << "\" is not present." << G4endl;
  return 0;
 }
 
 const RadmonSubPhysicsListLayout & subPhysicsList(subPhysicsListCollection.FindItemByLabel(physicsListName));
 return subPhysicsList.GetNAttributes();
}



const G4String &                                RadmonPhysicsLayout :: GetPhysicsListAttributeName(const G4String & physicsListName, G4int index) const
{
 if (!subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::GetPhysicsListAttributeName: Physics list \"" << physicsListName << "\" is not present." << G4endl;
  return GetNullStr();
 }
 
 const RadmonSubPhysicsListLayout & subPhysicsList(subPhysicsListCollection.FindItemByLabel(physicsListName));
 return subPhysicsList.GetAttributeName(index);
}



G4String                                        RadmonPhysicsLayout :: GetPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName, const G4String & defaultAttributeValue) const
{
 if (!subPhysicsListCollection.ExistsItemByLabel(physicsListName))
 {
  G4cout << "RadmonPhysicsLayout::GetPhysicsListAttributeName: Physics list \"" << physicsListName << "\" is not present." << G4endl;
  return GetNullStr();
 }
 
 const RadmonSubPhysicsListLayout & subPhysicsList(subPhysicsListCollection.FindItemByLabel(physicsListName));
 return subPhysicsList.GetAttribute(attributeName, defaultAttributeValue);
}





void                                            RadmonPhysicsLayout :: DumpLayout(std::ostream & out) const
{
 out << "- Physics list\n";
 const G4int n(subPhysicsListCollection.GetNItems());

 if (n==0)
 {
  out << "  - No physics list defined.\n";
  return;
 }

 const G4String indent("    - ");

 for(G4int i(0); i<n; i++)
 {
  out << "  - Physics list # " << i << '\n';
   
  subPhysicsListCollection.GetItem(i).DumpLayout(out, indent);
 }
}





G4bool                                          RadmonPhysicsLayout :: Load(std::istream & /* in */)
{
 // TO BE DONE
 G4cout << "RadmonPhysicsLayout::Load(): PLEASE CHECK" << G4endl;

 return false; 
}



G4bool                                          RadmonPhysicsLayout :: Save(std::ostream & /* out */) const
{
 // TO BE DONE
 G4cout << "RadmonPhysicsLayout::Save(): PLEASE CHECK" << G4endl;

 return false; 
}





inline G4String &                               RadmonPhysicsLayout :: GetNullStr() const
{
 static G4String null("");
 
 return null;
}



