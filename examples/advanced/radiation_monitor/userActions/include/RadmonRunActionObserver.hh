//
// File name:     RadmonRunActionObserver.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunActionObserver.hh,v 1.1 2005-11-24 02:31:56 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Run Action
//

#ifndef   RADMONRUNACTIONOBSERVER_HH
 #define  RADMONRUNACTIONOBSERVER_HH

 class G4Run;
 
 class RadmonRunActionObserver
 {
  public:
   inline                                       RadmonRunActionObserver();
   inline virtual                              ~RadmonRunActionObserver();
   
   virtual void                                 OnBeginOfRun(const G4Run * run) = 0;
   virtual void                                 OnEndOfRun(const G4Run * run) = 0;

  private:
                                                RadmonRunActionObserver(const RadmonRunActionObserver & copy);
   RadmonRunActionObserver &                    operator=(const RadmonRunActionObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonRunActionObserver.icc"
#endif /* RADMONRUNACTIONOBSERVER_HH */
