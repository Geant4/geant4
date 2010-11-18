/*
 * =============================================================================
 *
 *       Filename:  CexmcPhysicsManagerMessenger.hh
 *
 *    Description:  physics manager messenger (max IL correction etc.)
 *
 *        Version:  1.0
 *        Created:  16.10.2010 14:09:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PHYSICS_MANAGER_MESSENGER_HH
#define CEXMC_PHYSICS_MANAGER_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithADoubleAndUnit;
class  CexmcPhysicsManager;


class  CexmcPhysicsManagerMessenger : public G4UImessenger
{
    public:
        explicit CexmcPhysicsManagerMessenger(
                                        CexmcPhysicsManager *  physicsManager );

        ~CexmcPhysicsManagerMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcPhysicsManager *        physicsManager;

        G4UIcmdWithADoubleAndUnit *  setMaxILCorrection;
};


#endif

