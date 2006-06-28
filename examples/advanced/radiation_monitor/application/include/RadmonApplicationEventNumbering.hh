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
// File name:     RadmonApplicationEventNumbering.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventNumbering.hh,v 1.2 2006-06-28 13:45:06 gunter Exp $
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
