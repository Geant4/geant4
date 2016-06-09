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
 *       Filename:  CexmcPrimaryGeneratorActionMessenger.hh
 *
 *    Description:  user assigned gun parameters
 *
 *        Version:  1.0
 *        Created:  02.11.2009 13:08:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRIMARY_GENERATOR_ACTION_MESSENGER_HH
#define CEXMC_PRIMARY_GENERATOR_ACTION_MESSENGER_HH

#include <G4UImessenger.hh>

class  CexmcPrimaryGeneratorAction;
class  G4UIcommand;
class  G4UIcmdWithADouble;
class  G4UIcmdWithADoubleAndUnit;


class  CexmcPrimaryGeneratorActionMessenger : public G4UImessenger
{
    public:
        explicit CexmcPrimaryGeneratorActionMessenger(
                        CexmcPrimaryGeneratorAction *  primaryGeneratorAction );

        ~CexmcPrimaryGeneratorActionMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcPrimaryGeneratorAction *  primaryGeneratorAction;

        G4UIcmdWithADoubleAndUnit *    fwhmPosX;

        G4UIcmdWithADoubleAndUnit *    fwhmPosY;

        G4UIcmdWithADoubleAndUnit *    fwhmDirX;

        G4UIcmdWithADoubleAndUnit *    fwhmDirY;

        G4UIcmdWithADouble *           fwhmMomentumAmp;
};


#endif

