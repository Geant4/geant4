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
 *       Filename:  CexmcRunManagerMessenger.hh
 *
 *    Description:  init parameters (production model, gdml file etc.)
 *
 *        Version:  1.0
 *        Created:  03.11.2009 20:36:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RUN_MANAGER_MESSENGER_HH
#define CEXMC_RUN_MANAGER_MESSENGER_HH

#include <G4UImessenger.hh>

class  CexmcRunManager;
class  G4UIcommand;
class  G4UIcmdWithAString;
class  G4UIcmdWithAnInteger;
class  G4UIcmdWithABool;
class  G4UIcmdWithoutParameter;


class  CexmcRunManagerMessenger : public G4UImessenger
{
    public:
        explicit CexmcRunManagerMessenger( CexmcRunManager *  runManager );

        ~CexmcRunManagerMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcRunManager *          runManager;

        G4UIcmdWithAString *       setProductionModel;

        G4UIcmdWithAString *       setGdmlFile;

        G4UIcmdWithAString *       setGuiMacro;

        G4UIcmdWithAString *       setEventCountPolicy;

        G4UIcmdWithAString *       setEventDataVerboseLevel;

#ifdef CEXMC_USE_PERSISTENCY
        G4UIcmdWithAnInteger *     replayEvents;

        G4UIcmdWithAnInteger *     seekTo;

        G4UIcmdWithABool *         skipInteractionsWithoutEDT;
#endif

        G4UIcmdWithoutParameter *  registerScenePrimitives;

        G4UIcmdWithABool *         validateGdmlFile;
};


#endif

