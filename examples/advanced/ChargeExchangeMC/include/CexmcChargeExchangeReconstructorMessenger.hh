/*
 * =============================================================================
 *
 *       Filename:  CexmcChargeExchangeReconstructorMessenger.hh
 *
 *    Description:  charge exchange reconstructor messenger
 *
 *        Version:  1.0
 *        Created:  14.12.2009 17:49:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_CHARGE_EXCHANGE_RECONSTRUCTOR_MESSENGER_HH
#define CEXMC_CHARGE_EXCHANGE_RECONSTRUCTOR_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithABool;
class  G4UIcmdWithADoubleAndUnit;
class  G4String;
class  CexmcChargeExchangeReconstructor;


class  CexmcChargeExchangeReconstructorMessenger : public G4UImessenger
{
    public:
        explicit CexmcChargeExchangeReconstructorMessenger(
                            CexmcChargeExchangeReconstructor *  reconstructor );

        ~CexmcChargeExchangeReconstructorMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcChargeExchangeReconstructor *  reconstructor;

        G4UIcmdWithABool *                  useTableMass;

        G4UIcmdWithABool *                  useMassCut;

        G4UIcmdWithADoubleAndUnit *         mCutOPCenter;

        G4UIcmdWithADoubleAndUnit *         mCutNOPCenter;

        G4UIcmdWithADoubleAndUnit *         mCutOPWidth;

        G4UIcmdWithADoubleAndUnit *         mCutNOPWidth;

        G4UIcmdWithADoubleAndUnit *         mCutAngle;

        G4UIcmdWithABool *                  useAbsorbedEnergyCut;

        G4UIcmdWithADoubleAndUnit *         aeCutCLCenter;

        G4UIcmdWithADoubleAndUnit *         aeCutCRCenter;

        G4UIcmdWithADoubleAndUnit *         aeCutCLWidth;

        G4UIcmdWithADoubleAndUnit *         aeCutCRWidth;

        G4UIcmdWithADoubleAndUnit *         aeCutAngle;
};


#endif

