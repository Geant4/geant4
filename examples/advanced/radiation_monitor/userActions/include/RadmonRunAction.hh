//
// File name:     RadmonRunAction.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunAction.hh,v 1.1 2005-11-24 02:31:56 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon run action class
//

#ifndef   RADMONRUNACTION_HH
 #define  RADMONRUNACTION_HH
 
 // Include files
 #include "G4UserRunAction.hh"
 #include "globals.hh"
 #include <set>
 
 // Forward declarations
 class RadmonRunActionObserver;
 
 class RadmonRunAction : public G4UserRunAction
 {
  public:
   inline virtual                              ~RadmonRunAction();
   
   inline static RadmonRunAction *              Instance(void);

   void                                         AttachObserver(RadmonRunActionObserver * observer);
   void                                         DetachObserver(RadmonRunActionObserver * observer);
   
   virtual void                                 BeginOfRunAction(const G4Run * run);
   virtual void                                 EndOfRunAction(const G4Run * run);
   
  private:
                                                RadmonRunAction();
  // Hidden constructors and operators
                                                RadmonRunAction(const RadmonRunAction & copy);
   RadmonRunAction &                            operator=(const RadmonRunAction & copy);
   
  // Private attributes
   typedef std::set<RadmonRunActionObserver *>  ObserversSet;
   ObserversSet                                 observersSet;
   
   static RadmonRunAction *                     instance;
 };
 
 // Inline implementations
 #include "RadmonRunAction.icc"
#endif /* RADMONRUNACTION_HH */
