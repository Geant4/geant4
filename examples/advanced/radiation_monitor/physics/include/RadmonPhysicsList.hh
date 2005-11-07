//
// File name:     RadmonPhysicsList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsList.hh,v 1.1 2005-11-07 17:52:36 capra Exp $
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

 class RadmonPhysicsList : public G4VUserPhysicsList, public RadmonVLayoutObserver
 {
  public:
                                                RadmonPhysicsList(RadmonVPhysicsLayout * layout, RadmonVSubPhysicsListFactory * factory);
   virtual                                     ~RadmonPhysicsList();

   virtual void                                 OnLayoutChange(void);

  protected:
   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);

  private:
  // Private methods
   void                                         Destruct(void);

  // Hidden constructors and operators
                                                RadmonPhysicsList();
                                                RadmonPhysicsList(const RadmonPhysicsList & copy);
   RadmonPhysicsList &                          operator=(const RadmonPhysicsList & copy);

  // Private data types
   typedef std::list<RadmonVSubPhysicsList *>   SubPhysiscsLists;

  // Private attributes
   RadmonVPhysicsLayout *                       physicsLayout;
   RadmonVSubPhysicsListFactory *               constructorsFactory;
    
   SubPhysiscsLists                             subPhysiscsLists;
 };
#endif /* RADMONPHYSICSLIST_HH */
