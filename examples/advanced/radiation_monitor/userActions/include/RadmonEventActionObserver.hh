//
// File name:     RadmonEventActionObserver.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonEventActionObserver.hh,v 1.1 2005-11-24 02:31:56 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Event Action
//

#ifndef   RADMONEVENTACTIONOBSERVER_HH
 #define  RADMONEVENTACTIONOBSERVER_HH

 class G4Event;
 
 class RadmonEventActionObserver
 {
  public:
   inline                                       RadmonEventActionObserver();
   inline virtual                              ~RadmonEventActionObserver();
   
   virtual void                                 OnBeginOfEvent(const G4Event * event) = 0;
   virtual void                                 OnEndOfEvent(const G4Event * event) = 0;

  private:
                                                RadmonEventActionObserver(const RadmonEventActionObserver & copy);
   RadmonEventActionObserver &                  operator=(const RadmonEventActionObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonEventActionObserver.icc"
#endif /* RADMONEVENTACTIONOBSERVER_HH */
