//
// File name:     RadmonMessengersSupport.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMessengersSupport.hh,v 1.1 2005-09-19 19:39:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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
