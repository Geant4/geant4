/*
 * =============================================================================
 *
 *       Filename:  CexmcParticleGunMessenger.hh
 *
 *    Description:  original position and momentum of the incident beam particle
 *
 *        Version:  1.0
 *        Created:  15.12.2009 13:54:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PARTICLE_GUN_MESSENGER_HH
#define CEXMC_PARTICLE_GUN_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithAString;
class  G4UIcmdWithADoubleAndUnit;
class  G4UIcmdWith3Vector;
class  G4UIcmdWith3VectorAndUnit;
class  CexmcParticleGun;


class  CexmcParticleGunMessenger : public G4UImessenger
{
    public:
        explicit CexmcParticleGunMessenger( CexmcParticleGun *  particleGun );

        ~CexmcParticleGunMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cnd, G4String  value );

    private:
        CexmcParticleGun *           particleGun;

        G4UIcmdWithAString *         setParticle;

        G4UIcmdWith3VectorAndUnit *  setOrigPosition;

        G4UIcmdWith3Vector *         setOrigDirection;

        G4UIcmdWithADoubleAndUnit *  setOrigMomentumAmp;
};


#endif

