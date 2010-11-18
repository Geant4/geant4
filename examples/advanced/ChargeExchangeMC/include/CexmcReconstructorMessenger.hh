/*
 * =============================================================================
 *
 *       Filename:  CexmcReconstructorMessenger.hh
 *
 *    Description:  reconstructor messenger
 *
 *        Version:  1.0
 *        Created:  02.12.2009 15:33:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RECONSTRUCTOR_MESSENGER_HH
#define CEXMC_RECONSTRUCTOR_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithAString;
class  G4UIcmdWithABool;
class  G4UIcmdWithADoubleAndUnit;
class  G4String;
class  CexmcReconstructor;


class  CexmcReconstructorMessenger : public G4UImessenger
{
    public:
        explicit CexmcReconstructorMessenger( CexmcReconstructor *
                                                                reconstructor );

        ~CexmcReconstructorMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcReconstructor *   reconstructor;

        G4UIcmdWithAString *   setCalorimeterEntryPointDefinitionAlgorithm;

        G4UIcmdWithAString *   setCalorimeterEntryPointDepthDefinitionAlgorithm;

        G4UIcmdWithAString *   setCrystalSelectionAlgorithm;

        G4UIcmdWithABool *     useInnerRefCrystal;

        G4UIcmdWithADoubleAndUnit *  setCalorimeterEntryPointDepth;
};


#endif

