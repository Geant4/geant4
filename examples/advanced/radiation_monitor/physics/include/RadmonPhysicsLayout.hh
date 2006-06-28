//
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
//
// File name:     RadmonPhysicsLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsLayout.hh,v 1.3 2006-06-28 13:54:53 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class to keep track of the configured physics list
//

#ifndef   RADMONPHYSICSLAYOUT_HH
 #define  RADMONPHYSICSLAYOUT_HH
 
 // Include files
 #include "RadmonVPhysicsLayout.hh"

 #include "RadmonSubPhysicsListLayout.hh"
 #include "RadmonTLabelledCollection.hh"
 
 class RadmonPhysicsLayout : public RadmonVPhysicsLayout
 {
  public:
                                                RadmonPhysicsLayout();
                                               ~RadmonPhysicsLayout();

   virtual void                                 AddPhysicsList(const G4String & physicsListName);
   virtual void                                 RemovePhysicsList(const G4String & physicsListName);
   virtual G4int                                GetNPhysicsLists(void) const;
   virtual const G4String &                     GetPhysicsListName(G4int index) const;

   virtual void                                 SetPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName, const G4String & attributeValue);
   virtual void                                 ClearPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName);
   virtual G4int                                GetPhysicsListNAttributes(const G4String & physicsListName) const;
   virtual const G4String &                     GetPhysicsListAttributeName(const G4String & physicsListName, G4int index) const;
   virtual G4String                             GetPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName, const G4String & defaultAttributeValue=G4String()) const;

   virtual void                                 DumpLayout(std::ostream & out) const;

   virtual G4bool                               Load(std::istream & in);
   virtual G4bool                               Save(std::ostream & out) const;

  private:
   inline G4String &                            GetNullStr() const;

  // Hidden constructors and operators
                                                RadmonPhysicsLayout(const RadmonPhysicsLayout & copy);
   RadmonPhysicsLayout &                        operator=(const RadmonPhysicsLayout & copy);

   RadmonTLabelledCollection<RadmonSubPhysicsListLayout> subPhysicsListCollection;
 };
#endif /* RADMONPHYSICSLAYOUT_HH */
