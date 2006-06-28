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
// File name:     RadmonVSubPhysicsList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsList.hh,v 1.3 2006-06-28 13:55:53 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
