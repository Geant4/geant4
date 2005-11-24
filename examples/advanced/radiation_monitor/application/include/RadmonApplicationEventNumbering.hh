//
// File name:     RadmonApplicationEventNumbering.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventNumbering.hh,v 1.1 2005-11-24 02:34:21 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Dumps event id
//

#ifndef   RADMONAPPLICATIONEVENTNUMBERING_HH
 #define  RADMONAPPLICATIONEVENTNUMBERING_HH
 
 // Include files
 #include "RadmonEventActionObserver.hh"
 #include "globals.hh"
 
 class RadmonApplicationEventNumbering : public RadmonEventActionObserver
 {
  public:
   inline                                       RadmonApplicationEventNumbering();
   inline virtual                              ~RadmonApplicationEventNumbering();
   
   inline void                                  SetDumpEvery(G4int events);
   
   virtual void                                 OnBeginOfEvent(const G4Event * event);
   inline virtual void                          OnEndOfEvent(const G4Event * event);
   
  private:
                                                RadmonApplicationEventNumbering(const RadmonApplicationEventNumbering & copy);
   RadmonApplicationEventNumbering &            operator=(const RadmonApplicationEventNumbering & copy);
   
   G4int                                        dumpEvery;
 };
 
 // Inline implementations
 #include "RadmonApplicationEventNumbering.icc"
#endif /* RADMONAPPLICATIONEVENTNUMBERING_HH */
