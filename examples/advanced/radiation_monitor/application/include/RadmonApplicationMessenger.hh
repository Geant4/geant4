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
// File name:     RadmonApplicationMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationMessenger.hh,v 1.1.2.2 2006/06/29 16:08:15 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   UI commands for managing application level options
//

#ifndef   RADMONAPPLICATIONMESSENGER_HH
 #define  RADMONAPPLICATIONMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonApplicationEventNumbering;
 class RadmonApplicationEventTracks;
 class RadmonApplicationRunNumbering;
 
 class RadmonApplicationMessenger : public RadmonMessenger
 {
  public:
                                                RadmonApplicationMessenger();
                                               ~RadmonApplicationMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonApplicationMessenger(const RadmonApplicationMessenger & copy);
   RadmonApplicationMessenger &                 operator=(const RadmonApplicationMessenger & copy);

  // Commands
   RadmonApplicationEventNumbering *            eventNumbering;
   RadmonApplicationEventTracks *               eventTracks;
   RadmonApplicationRunNumbering *              runNumbering;

   RADMON_DECLARE_COMMAND(EnableRunsDump);
   RADMON_DECLARE_COMMAND(DisableRunsDump);
   RADMON_DECLARE_COMMAND(DumpEventsEvery);
   RADMON_DECLARE_COMMAND(DisableEventsDump);
   RADMON_DECLARE_COMMAND(EnableTracksVisualisation);
   RADMON_DECLARE_COMMAND(DisableTracksVisualisation);
 };
#endif /* RADMONAPPLICATIONMESSENGER_HH */
