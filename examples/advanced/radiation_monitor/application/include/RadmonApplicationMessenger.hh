//
// File name:     RadmonApplicationMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationMessenger.hh,v 1.2 2006-01-06 12:52:31 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
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
   RADMON_DECLARE_COMMAND(SetSeed);
 };
#endif /* RADMONAPPLICATIONMESSENGER_HH */
