//
// File name:     RadmonSteppingActionObserver.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSteppingActionObserver.hh,v 1.1 2005-12-06 19:36:23 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Stepping Action
//

#ifndef   RADMONSTEPPINGACTIONOBSERVER_HH
 #define  RADMONSTEPPINGACTIONOBSERVER_HH

 class G4Step;
 
 class RadmonSteppingActionObserver
 {
  public:
   inline                                       RadmonSteppingActionObserver();
   inline virtual                              ~RadmonSteppingActionObserver();
   
   virtual void                                 OnUserStepping(const G4Step * step) = 0;

  private:
                                                RadmonSteppingActionObserver(const RadmonSteppingActionObserver & copy);
   RadmonSteppingActionObserver &               operator=(const RadmonSteppingActionObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonSteppingActionObserver.icc"
#endif /* RADMONSTEPPINGACTIONOBSERVER_HH */
