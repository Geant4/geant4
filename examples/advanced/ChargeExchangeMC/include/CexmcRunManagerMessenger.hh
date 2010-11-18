/*
 * =============================================================================
 *
 *       Filename:  CexmcRunManagerMessenger.hh
 *
 *    Description:  init parameters (production model, gdml file etc.)
 *
 *        Version:  1.0
 *        Created:  03.11.2009 20:36:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RUN_MANAGER_MESSENGER_HH
#define CEXMC_RUN_MANAGER_MESSENGER_HH

#include <G4UImessenger.hh>

class  CexmcRunManager;
class  G4UIcommand;
class  G4UIcmdWithAString;
class  G4UIcmdWithAnInteger;
class  G4UIcmdWithABool;


class  CexmcRunManagerMessenger : public G4UImessenger
{
    public:
        explicit CexmcRunManagerMessenger( CexmcRunManager *  runManager );

        ~CexmcRunManagerMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcRunManager *       runManager;

        G4UIcmdWithAString *    setProductionModel;

        G4UIcmdWithAString *    setGdmlFile;

        G4UIcmdWithAString *    setGuiMacro;

        G4UIcmdWithAString *    setEventCountPolicy;

        G4UIcmdWithAString *    setEventDataVerboseLevel;

        G4UIcmdWithAnInteger *  replayEvents;

        G4UIcmdWithAnInteger *  seekTo;

        G4UIcmdWithABool *      skipInteractionsWithoutEDT;

        G4UIcmdWithABool *      validateGdmlFile;
};


#endif

