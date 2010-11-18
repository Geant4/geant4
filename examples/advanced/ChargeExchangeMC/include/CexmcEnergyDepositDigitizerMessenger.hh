/*
 * =============================================================================
 *
 *       Filename:  CexmcEnergyDepositDigitizerMessenger.hh
 *
 *    Description:  energy deposit digitizer messenger
 *
 *        Version:  1.0
 *        Created:  29.11.2009 18:54:33
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_ENERGY_DEPOSIT_DIGITIZER_MESSENGER_HH
#define CEXMC_ENERGY_DEPOSIT_DIGITIZER_MESSENGER_HH

#include <G4UImessenger.hh>

class  CexmcEnergyDepositDigitizer;
class  G4UIcommand;
class  G4UIcmdWithADouble;
class  G4UIcmdWithADoubleAndUnit;
class  G4UIcmdWithAString;
class  G4UIcmdWithABool;
class  G4UIcmdWith3Vector;
class  G4UIcmdWithoutParameter;


class  CexmcEnergyDepositDigitizerMessenger : public G4UImessenger
{
    public:
        explicit CexmcEnergyDepositDigitizerMessenger(
                        CexmcEnergyDepositDigitizer *  energyDepositDigitiser );

        ~CexmcEnergyDepositDigitizerMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcEnergyDepositDigitizer *  energyDepositDigitizer;

        G4UIcmdWithADoubleAndUnit *    setMonitorThreshold;

        G4UIcmdWithADoubleAndUnit *    setVetoCountersThreshold;

        G4UIcmdWithADoubleAndUnit *    setLeftVetoCounterThreshold;

        G4UIcmdWithADoubleAndUnit *    setRightVetoCounterThreshold;

        G4UIcmdWithADoubleAndUnit *    setCalorimetersThreshold;

        G4UIcmdWithADoubleAndUnit *    setLeftCalorimeterThreshold;

        G4UIcmdWithADoubleAndUnit *    setRightCalorimeterThreshold;

        G4UIcmdWithAString *           setCalorimeterTriggerAlgorithm;

        G4UIcmdWithAString *           setOuterCrystalsVetoAlgorithm;

        G4UIcmdWithADouble *           setOuterCrystalsVetoFraction;

        G4UIcmdWithABool *             applyFiniteCrystalResolution;

        G4UIcmdWith3Vector *           addCrystalResolutionRange;

        G4UIcmdWithoutParameter *      clearCrystalResolutionData;
};


#endif

