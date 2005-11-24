//
// File name:     RadmonApplicationEventTracks.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventTracks.hh,v 1.1 2005-11-24 02:34:21 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Dumps event id
//

#ifndef   RADMONAPPLICATIONEVENTTRACKS_HH
 #define  RADMONAPPLICATIONEVENTTRACKS_HH
 
 // Include files
 #include "RadmonEventActionObserver.hh"
 #include "globals.hh"
 
 class RadmonApplicationEventTracks : public RadmonEventActionObserver
 {
  public:
   inline                                       RadmonApplicationEventTracks();
   inline virtual                              ~RadmonApplicationEventTracks();
   
   inline void                                  Enable(void);
   inline void                                  Disable(void);
   
   inline virtual void                          OnBeginOfEvent(const G4Event * event);
   virtual void                                 OnEndOfEvent(const G4Event * event);
   
  private:
                                                RadmonApplicationEventTracks(const RadmonApplicationEventTracks & copy);
   RadmonApplicationEventTracks &            operator=(const RadmonApplicationEventTracks & copy);

   G4bool                                       enable;
 };
 
 // Inline implementations
 #include "RadmonApplicationEventTracks.icc"
#endif /* RADMONAPPLICATIONEVENTTRACKS_HH */
