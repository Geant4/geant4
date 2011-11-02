//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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
class  G4UIcmdWithAString;
class  G4UIcmdWithADoubleAndUnit;
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

        G4UIcmdWithADoubleAndUnit *         setExpectedMomentumAmp;

        G4UIcmdWithADoubleAndUnit *         setExpectedMomentumAmpDiff;

        G4UIcmdWithAString *                setEDCollectionAlgorithm;
};


#endif

