//
// File name:     RadmonEventAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonEventAction.hh,v 1.1 2005-11-24 02:31:56 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon event action class
//

#ifndef   RADMONEVENTACTION_HH
 #define  RADMONEVENTACTION_HH
 
 // Include files
 #include "G4UserEventAction.hh"
 #include "globals.hh"
 #include <set>
 
 // Forward declarations
 class RadmonEventActionObserver;
 
 class RadmonEventAction : public G4UserEventAction
 {
  public:
   inline virtual                              ~RadmonEventAction();
   
   inline static RadmonEventAction *            Instance(void);

   void                                         AttachObserver(RadmonEventActionObserver * observer);
   void                                         DetachObserver(RadmonEventActionObserver * observer);
   
   virtual void                                 BeginOfEventAction(const G4Event * event);
   virtual void                                 EndOfEventAction(const G4Event * event);
   
  private:
                                                RadmonEventAction();
  // Hidden constructors and operators
                                                RadmonEventAction(const RadmonEventAction & copy);
   RadmonEventAction &                          operator=(const RadmonEventAction & copy);
   
  // Private attributes
   typedef std::set<RadmonEventActionObserver *> ObserversSet;
   ObserversSet                                 observersSet;
   
   static RadmonEventAction *                   instance;
 };
 
 // Inline implementations
 #include "RadmonEventAction.icc"
#endif /* RADMONEVENTACTION_HH */
