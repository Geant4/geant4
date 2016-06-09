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
//
// File name:     RadmonMessengersSupport.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMessengersSupport.hh,v 1.3 2006/06/29 16:14:27 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Define macros useful to maintain the code of messengers
//

#ifndef   RADMONMESSENGERSSUPPORT_HH
 #define  RADMONMESSENGERSSUPPORT_HH

 #include "G4UIcommand.hh"

 #define RADMON_DECLARE_COMMAND(command)        private:                                      \
                                                 void On ## command (const G4String & value); \
                                                 G4UIcommand * cmd ## command
 
 #define RADMON_INITIALIZE_COMMAND(command)     cmd ## command (0)

 #define _CREATE_COMMAND(command, guidance)     cmd ## command = new G4UIcommand(COMMANDS_PATH #command , this);                                                 \
                                                                                                                                                                 \
                                                if (! cmd ## command )                                                                                           \
                                                {                                                                                                                \
                                                 G4cout << "RadmonDetectorMessenger::RadmonDetectorMessenger: Command \"" #command "\" not allocated."<< G4endl; \
                                                 return;                                                                                                         \
                                                }                                                                                                                \
                                                                                                                                                                 \
                                                cmd ## command ->SetGuidance(guidance)

 #define _ADD_ARG(command, argName)              cmd ## command ->SetParameter(new G4UIparameter(argName, 's', false))

 #define RADMON_CREATE_COMMAND_0ARGS(command, guidance)                             \
                                                _CREATE_COMMAND(command, guidance);

 #define RADMON_CREATE_COMMAND_1ARG(command, guidance, argName0)                                \
                                                RADMON_CREATE_COMMAND_0ARGS(command, guidance); \
                                                _ADD_ARG(command, argName0)

 #define RADMON_CREATE_COMMAND_2ARGS(command, guidance, argName0, argName1)                              \
                                                RADMON_CREATE_COMMAND_1ARG(command, guidance, argName0); \
                                                _ADD_ARG(command, argName1)

 #define RADMON_CREATE_COMMAND_3ARGS(command, guidance, argName0, argName1, argName2)                               \
                                                RADMON_CREATE_COMMAND_2ARGS(command, guidance, argName0, argName1); \
                                                _ADD_ARG(command, argName2)

 #define RADMON_CREATE_COMMAND_4ARGS(command, guidance, argName0, argName1, argName2, argName3)                               \
                                                RADMON_CREATE_COMMAND_3ARGS(command, guidance, argName0, argName1, argName2); \
                                                _ADD_ARG(command, argName3)

 #define RADMON_CREATE_COMMAND_5ARGS(command, guidance, argName0, argName1, argName2, argName3, argName4)                               \
                                                RADMON_CREATE_COMMAND_4ARGS(command, guidance, argName0, argName1, argName2, argName3); \
                                                _ADD_ARG(command, argName4)

 #define RADMON_CREATE_COMMAND_6ARGS(command, guidance, argName0, argName1, argName2, argName3, argName4, argName5)                               \
                                                RADMON_CREATE_COMMAND_5ARGS(command, guidance, argName0, argName1, argName2, argName3, argName4); \
                                                _ADD_ARG(command, argName5)

 #define RADMON_DESTROY_COMMAND(command)        if ( cmd ## command )   \
                                                 delete cmd ## command;
                                                
 #define RADMON_BEGIN_LIST_SET_COMMANDS         if (!command)                                                             \
                                                 G4cout << "RadmonDetectorMessenger::SetNewValue: command==0." << G4endl; \
                                                else

 #define  RADMON_SET_COMMAND(name)              if ( cmd ## name == command) \
                                                 On ## name (newValue);      \
                                                else

 #define RADMON_END_LIST_SET_COMMANDS           G4cout << "RadmonDetectorMessenger::SetNewValue: Command \"" << command->GetCommandPath() << "\" not supported." << G4endl;
#endif /* RADMONMESSENGERSSUPPORT_HH */
