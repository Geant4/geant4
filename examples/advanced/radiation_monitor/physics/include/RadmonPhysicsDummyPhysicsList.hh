//
// File name:     RadmonPhysicsDummyPhysicsList.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsDummyPhysicsList.hh,v 1.1 2005-09-14 12:31:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Dummy physics list
//

#ifndef   RADMONPHYSICSDUMMYPHYSICSLIST_HH
 #define  RADMONPHYSICSDUMMYPHYSICSLIST_HH
 
 //Include files
 #include "G4VUserPhysicsList.hh"

 class RadmonPhysicsDummyPhysicsList : public G4VUserPhysicsList
 {
  public:
   inline                                       RadmonPhysicsDummyPhysicsList();
   inline virtual                              ~RadmonPhysicsDummyPhysicsList();

  protected:
   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);

  private:
  // Hidden constructors and operators
                                                RadmonPhysicsDummyPhysicsList(const RadmonPhysicsDummyPhysicsList & copy);
   RadmonPhysicsDummyPhysicsList &              operator=(const RadmonPhysicsDummyPhysicsList & copy);
 };
 
 // Inline implementations
 #include "RadmonPhysicsDummyPhysicsList.icc"
#endif /* RADMONPHYSICSDUMMYPHYSICSLIST_HH */
