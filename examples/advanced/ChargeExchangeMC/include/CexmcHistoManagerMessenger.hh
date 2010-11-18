/*
 * =============================================================================
 *
 *       Filename:  CexmcHistoManagerMessenger.hh
 *
 *    Description:  commands to list and show histograms
 *
 *        Version:  1.0
 *        Created:  17.12.2009 21:38:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_HISTO_MANAGER_MESSENGER_HH
#define CEXMC_HISTO_MANAGER_MESSENGER_HH

#ifdef CEXMC_USE_ROOT

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithoutParameter;
class  G4UIcmdWithAString;


class  CexmcHistoManagerMessenger : public G4UImessenger
{
    public:
        CexmcHistoManagerMessenger();

        ~CexmcHistoManagerMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        G4UIcmdWithoutParameter *  listHistos;

        G4UIcmdWithAString *       printHisto;

#ifdef CEXMC_USE_ROOTQT
        G4UIcmdWithAString *       drawHisto;
#endif
};

#endif


#endif

