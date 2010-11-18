/*
 * =============================================================================
 *
 *       Filename:  CexmcPrimaryGeneratorActionMessenger.hh
 *
 *    Description:  user assigned gun parameters
 *
 *        Version:  1.0
 *        Created:  02.11.2009 13:08:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRIMARY_GENERATOR_ACTION_MESSENGER_HH
#define CEXMC_PRIMARY_GENERATOR_ACTION_MESSENGER_HH

#include <G4UImessenger.hh>

class  CexmcPrimaryGeneratorAction;
class  G4UIcommand;
class  G4UIcmdWithADouble;
class  G4UIcmdWithADoubleAndUnit;


class  CexmcPrimaryGeneratorActionMessenger : public G4UImessenger
{
    public:
        explicit CexmcPrimaryGeneratorActionMessenger(
                        CexmcPrimaryGeneratorAction *  primaryGeneratorAction );

        ~CexmcPrimaryGeneratorActionMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcPrimaryGeneratorAction *  primaryGeneratorAction;

        G4UIcmdWithADoubleAndUnit *    fwhmPosX;

        G4UIcmdWithADoubleAndUnit *    fwhmPosY;

        G4UIcmdWithADoubleAndUnit *    fwhmDirX;

        G4UIcmdWithADoubleAndUnit *    fwhmDirY;

        G4UIcmdWithADouble *           fwhmMomentumAmp;
};


#endif

