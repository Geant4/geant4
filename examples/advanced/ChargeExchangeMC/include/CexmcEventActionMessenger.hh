/*
 * =============================================================================
 *
 *       Filename:  CexmcEventActionMessenger.hh
 *
 *    Description:  event action messenger (verbose level etc.)
 *
 *        Version:  1.0
 *        Created:  25.11.2009 14:38:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_EVENT_ACTION_MESSENGER_HH
#define CEXMC_EVENT_ACTION_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithAnInteger;
class  G4UIcmdWithABool;
class  G4String;
class  CexmcEventAction;


class  CexmcEventActionMessenger : public G4UImessenger
{
    public:
        explicit CexmcEventActionMessenger( CexmcEventAction *  eventAction );

        ~CexmcEventActionMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcEventAction *      eventAction;

        G4UIcmdWithAnInteger *  setVerboseLevel;

        G4UIcmdWithAnInteger *  setVerboseDrawLevel;

        G4UIcmdWithABool *      drawTrajectoryMarkers;
};


#endif

