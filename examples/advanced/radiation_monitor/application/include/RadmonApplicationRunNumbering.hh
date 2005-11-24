//
// File name:     RadmonApplicationRunNumbering.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationRunNumbering.hh,v 1.1 2005-11-24 02:34:21 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Dumps run id
//

#ifndef   RADMONAPPLICATIONRUNNUMBERING_HH
 #define  RADMONAPPLICATIONRUNNUMBERING_HH
 
 // Include files
 #include "RadmonRunActionObserver.hh"
 #include "globals.hh"
 
 class RadmonApplicationRunNumbering : public RadmonRunActionObserver
 {
  public:
   inline                                       RadmonApplicationRunNumbering();
   inline virtual                              ~RadmonApplicationRunNumbering();
   
   inline void                                  Enable(void);
   inline void                                  Disable(void);
   
   virtual void                                 OnBeginOfRun(const G4Run * run);
   inline virtual void                          OnEndOfRun(const G4Run * run);
   
  private:
                                                RadmonApplicationRunNumbering(const RadmonApplicationRunNumbering & copy);
   RadmonApplicationRunNumbering &              operator=(const RadmonApplicationRunNumbering & copy);
   
   G4bool                                       enable;
 };
 
 // Inline implementations
 #include "RadmonApplicationRunNumbering.icc"
#endif /* RADMONAPPLICATIONRUNNUMBERING_HH */
