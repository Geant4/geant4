//
// File name:     RadmonPhysicsSteppingAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsSteppingAction.hh,v 1.1 2005-12-06 19:36:23 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Stepping Action
//

#ifndef   RADMONPHYSICSSTEPPINGACTION_HH
 #define  RADMONPHYSICSSTEPPINGACTION_HH

 #include "RadmonSteppingActionObserver.hh"
 
 class RadmonPhysicsSteppingAction : public RadmonSteppingActionObserver
 {
  public:
   inline                                       RadmonPhysicsSteppingAction();
   inline virtual                              ~RadmonPhysicsSteppingAction();
   
   virtual void                                 OnUserStepping(const G4Step * step);

  private:
                                                RadmonPhysicsSteppingAction(const RadmonPhysicsSteppingAction & copy);
   RadmonPhysicsSteppingAction &                operator=(const RadmonPhysicsSteppingAction & copy);
 };
 
 // Inline implementations
 #include "RadmonPhysicsSteppingAction.icc"
#endif /* RADMONPHYSICSSTEPPINGACTION_HH */
