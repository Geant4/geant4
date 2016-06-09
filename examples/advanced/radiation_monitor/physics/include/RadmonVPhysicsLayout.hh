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
// File name:     RadmonVPhysicsLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVPhysicsLayout.hh,v 1.4 2006/06/29 16:18:23 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
