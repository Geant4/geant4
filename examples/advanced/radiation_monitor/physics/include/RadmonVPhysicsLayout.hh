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
// File name:     RadmonVPhysicsLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVPhysicsLayout.hh,v 1.3 2006-06-28 13:55:49 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class to keep track of the configured physics list
//

#ifndef   RADMONVPHYSICSLAYOUT_HH
 #define  RADMONVPHYSICSLAYOUT_HH
 
 // Include files
 #include "RadmonVLayoutSubject.hh"

 #include "globals.hh"
 
 class RadmonVPhysicsLayout : public RadmonVLayoutSubject
 {
  public:
   virtual void                                 AddPhysicsList(const G4String & physicsListName) = 0;
   virtual void                                 RemovePhysicsList(const G4String & physicsListName) = 0;
   virtual G4int                                GetNPhysicsLists(void) const = 0;
   virtual const G4String &                     GetPhysicsListName(G4int index) const = 0;

   virtual void                                 SetPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName, const G4String & attributeValue) = 0;
   virtual void                                 ClearPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName) = 0;
   virtual G4int                                GetPhysicsListNAttributes(const G4String & physicsListName) const = 0;
   virtual const G4String &                     GetPhysicsListAttributeName(const G4String & physicsListName, G4int index) const = 0;
   virtual G4String                             GetPhysicsListAttribute(const G4String & physicsListName, const G4String & attributeName, const G4String & defaultAttributeValue=G4String()) const = 0;

   virtual void                                 DumpLayout(std::ostream & out) const = 0;

   virtual G4bool                               Load(std::istream & in) = 0;
   virtual G4bool                               Save(std::ostream & out) const = 0;

  protected:
   inline                                       RadmonVPhysicsLayout();
   inline                                      ~RadmonVPhysicsLayout();

  private:
  // Hidden constructors and operators
                                                RadmonVPhysicsLayout(const RadmonVPhysicsLayout & copy);
    RadmonVPhysicsLayout &                      operator=(const RadmonVPhysicsLayout & copy);
 };
 
 // Inline implementations
 #include "RadmonVPhysicsLayout.icc"
#endif /* RADMONVPHYSICSLAYOUT_HH */
