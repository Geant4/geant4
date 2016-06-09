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
// File name:     RadmonVSubPhysicsList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsList.hh,v 1.4 2006/06/29 16:18:27 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Abstract class of a physics list piece
//

#ifndef   RADMONVSUBPHYSICSLIST_HH
 #define  RADMONVSUBPHYSICSLIST_HH
 
 // Forward declarations
 class G4LogicalVolume;
 class RadmonPhysicsInfoList;
 class G4String;
 
 class RadmonVSubPhysicsList
 {
  public:
   inline virtual                              ~RadmonVSubPhysicsList();
    
   virtual void                                 SetPhysicsListAttribute(const G4String & attributeName, const G4String &value) = 0;

   virtual void                                 ConstructParticle(void) = 0;
   virtual void                                 ConstructProcess(void) = 0;
   virtual void                                 SetCuts(void) = 0;
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const = 0 ;

  protected:
   inline                                       RadmonVSubPhysicsList();

  private:
  // Hidden constructors and operators
                                                RadmonVSubPhysicsList(const RadmonVSubPhysicsList & copy);
   RadmonVSubPhysicsList &                      operator=(const RadmonVSubPhysicsList & copy);
 };
 
 // Inline implementations
 #include "RadmonVSubPhysicsList.icc"
#endif /* RADMONVSUBPHYSICSLIST_HH */
