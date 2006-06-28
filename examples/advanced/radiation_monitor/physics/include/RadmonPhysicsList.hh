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
// File name:     RadmonPhysicsList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsList.hh,v 1.5 2006-06-28 13:54:55 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Implementation of the G4VUserPhysicsList
//

#ifndef   RADMONPHYSICSLIST_HH
 #define  RADMONPHYSICSLIST_HH

 // Include files
 #include "globals.hh"
 #include "G4VUserPhysicsList.hh"
 #include "RadmonVLayoutObserver.hh"
 #include <list>
 
 // Forward declaration
 class RadmonVPhysicsLayout;
 class RadmonVSubPhysicsListFactory;
 class RadmonVSubPhysicsList;
 class RadmonPhysicsSteppingAction;

 class RadmonPhysicsList : public G4VUserPhysicsList, public RadmonVLayoutObserver
 {
  public:
                                                RadmonPhysicsList(RadmonVPhysicsLayout * layout, RadmonVSubPhysicsListFactory * factory);
   virtual                                     ~RadmonPhysicsList();

   virtual void                                 OnLayoutChange(void);

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);

  private:
  // Private methods
   void                                         Destruct(void);
   void                                         CheckUpdate(void);
   void                                         UpdateProcessManagers(void);

  // Hidden constructors and operators
                                                RadmonPhysicsList();
                                                RadmonPhysicsList(const RadmonPhysicsList & copy);
   RadmonPhysicsList &                          operator=(const RadmonPhysicsList & copy);

  // Private data types
   typedef std::list<RadmonVSubPhysicsList *>   SubPhysiscsLists;

  // Private attributes
   RadmonVPhysicsLayout *                       physicsLayout;
   RadmonVSubPhysicsListFactory *               subPhysicsListFactory;
    
   SubPhysiscsLists                             subPhysiscsLists;

   RadmonPhysicsSteppingAction *                steppingAction;
   
   G4bool                                       initializationMethodsCalled;
   G4bool                                       changed;
 };
#endif /* RADMONPHYSICSLIST_HH */
