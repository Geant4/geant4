/*
 * ============================================================================
 *
 *       Filename:  CexmcSteppingAction.hh
 *
 *    Description:  stepping action
 *
 *        Version:  1.0
 *        Created:  27.10.2009 15:59:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifndef CEXMC_STEPPING_ACTION_HH
#define CEXMC_STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>

class  G4Step;
class  G4LogicalVolume;
class  CexmcPhysicsManager;


class  CexmcSteppingAction : public G4UserSteppingAction
{
    public:
        explicit CexmcSteppingAction( CexmcPhysicsManager *  physicsManager );

    public:
        void  UserSteppingAction( const G4Step *  step );

    private:
        CexmcPhysicsManager *    physicsManager;

        const G4LogicalVolume *  targetVolume;
};


#endif

