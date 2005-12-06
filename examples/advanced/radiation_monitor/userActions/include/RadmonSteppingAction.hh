//
// File name:     RadmonSteppingAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingAction.hh,v 1.1 2005-12-06 19:36:23 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon stepping action class
//

#ifndef   RADMONSTEPPINGACTION_HH
 #define  RADMONSTEPPINGACTION_HH
 
 // Include files
 #include "G4UserSteppingAction.hh"
 #include "globals.hh"
 #include <set>
 
 // Forward declarations
 class RadmonSteppingActionObserver;
 
 class RadmonSteppingAction : public G4UserSteppingAction
 {
  public:
   inline virtual                              ~RadmonSteppingAction();
   
   inline static RadmonSteppingAction *         Instance(void);

   void                                         AttachObserver(RadmonSteppingActionObserver * observer);
   void                                         DetachObserver(RadmonSteppingActionObserver * observer);
   
   virtual void                                 UserSteppingAction(const G4Step * step);
   
  private:
                                                RadmonSteppingAction();
  // Hidden constructors and operators
                                                RadmonSteppingAction(const RadmonSteppingAction & copy);
   RadmonSteppingAction &                       operator=(const RadmonSteppingAction & copy);
   
  // Private attributes
   typedef std::set<RadmonSteppingActionObserver *> ObserversSet;
   ObserversSet                                 observersSet;
   
   static RadmonSteppingAction *                instance;
 };
 
 // Inline implementations
 #include "RadmonSteppingAction.icc"
#endif /* RADMONSTEPPINGACTION_HH */
