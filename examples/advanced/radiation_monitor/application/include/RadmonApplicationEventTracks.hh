//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonApplicationEventTracks.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventTracks.hh,v 1.3 2006/06/29 16:08:07 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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
