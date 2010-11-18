/*
 * =============================================================================
 *
 *       Filename:  CexmcProductionModelMessenger.hh
 *
 *    Description:  set various production model aspects
 *
 *        Version:  1.0
 *        Created:  03.11.2009 15:50:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRODUCTION_MODEL_MESSENGER_HH
#define CEXMC_PRODUCTION_MODEL_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithABool;
class  G4UIcmdWith3Vector;
class  CexmcProductionModel;


class  CexmcProductionModelMessenger : public G4UImessenger
{
    public:
        explicit CexmcProductionModelMessenger( CexmcProductionModel *
                                                productionModel );

        ~CexmcProductionModelMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcProductionModel *  productionModel;

        G4UIcmdWithABool *      applyFermiMotion;

        G4UIcmdWith3Vector *    setAngularRange;

        G4UIcmdWith3Vector *    addAngularRange;
};


#endif

