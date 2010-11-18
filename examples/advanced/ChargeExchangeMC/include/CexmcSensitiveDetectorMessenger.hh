/*
 * =============================================================================
 *
 *       Filename:  CexmcSensitiveDetectorMessenger.hh
 *
 *    Description:  sensitive detector messenger (verbose level etc.)
 *
 *        Version:  1.0
 *        Created:  15.11.2009 14:03:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SENSITIVE_DETECTOR_MESSENGER_HH
#define CEXMC_SENSITIVE_DETECTOR_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4VPrimitiveScorer;
class  G4UIcommand;
class  G4UIcmdWithAnInteger;
class  G4UIdirectory;
class  G4String;


class  CexmcSensitiveDetectorMessenger : public G4UImessenger
{
    public:
        CexmcSensitiveDetectorMessenger( G4VPrimitiveScorer *  scorer,
                                         const G4String &  detectorName );

        ~CexmcSensitiveDetectorMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        G4VPrimitiveScorer *    scorer;

        G4UIdirectory *         detectorPath;

        G4UIcmdWithAnInteger *  setVerboseLevel;
};


#endif

