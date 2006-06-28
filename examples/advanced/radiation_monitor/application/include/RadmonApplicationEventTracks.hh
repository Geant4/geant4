//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonApplicationEventTracks.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventTracks.hh,v 1.2 2006-06-28 13:45:10 gunter Exp $
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
