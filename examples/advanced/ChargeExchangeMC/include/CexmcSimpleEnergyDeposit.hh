/*
 * =============================================================================
 *
 *       Filename:  CexmcSimpleEnergyDeposit.hh
 *
 *    Description:  simple energy deposit scorer
 *
 *        Version:  1.0
 *        Created:  14.11.2009 12:45:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_ENERGY_DEPOSIT_HH
#define CEXMC_SIMPLE_ENERGY_DEPOSIT_HH

#include <G4VPrimitiveScorer.hh>
#include <G4THitsMap.hh>

class  G4HCofThisEvent;
class  G4Step;
class  CexmcSensitiveDetectorMessenger;


typedef G4THitsMap< G4double >  CexmcEnergyDepositCollection;


class  CexmcSimpleEnergyDeposit : public G4VPrimitiveScorer
{
    public:
        explicit CexmcSimpleEnergyDeposit( const G4String &  name );

        virtual ~CexmcSimpleEnergyDeposit();

    public:
        void   Initialize( G4HCofThisEvent *  hcOfThisEvent );

        void   EndOfEvent( G4HCofThisEvent *  hcOfThisEvent );

        void   DrawAll( void );

        void   PrintAll( void );

        void   clear( void );

    protected:
        G4int   GetIndex( G4Step *  step );

        G4bool  ProcessHits( G4Step *  step, G4TouchableHistory *  tHistory );

    protected:
        CexmcEnergyDepositCollection *     eventMap;

    private:
        CexmcSensitiveDetectorMessenger *  messenger;

        G4int                              hcId;
};


#endif

