// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIxvt.cc,v 1.2 1999-04-13 01:25:02 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifdef G4UI_BUILD_XVT_SESSION

////////////////////////////////////////////////////////////////////////////
//                   XVT driver class implementation
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Written by: Simon Prior
//       Date: 23/04/97
//
// Updated for state machine: 12/08/97
//
// The update consisted of adding the new state methods and I also deleted
// a lot of old redundant code (like ChangeDirectory etc which aren't 
// needed with the XVT GUI).
//
////////////////////////////////////////////////////////////////////////////
#include "G4UIxvt.hh"
#include "G4StateManager.hh"
#include "G4UIxvtMessenger.hh"

////////////////////////////////////////////////////////////////////////////
//                   
// G4UIxvt - Consructor
//
// Initalise the class variables.
//
// Note: number is the number of bytes transfered with each read or write.
//
////////////////////////////////////////////////////////////////////////////
G4UIxvt::G4UIxvt()
{  
  UI = G4UImanager::GetUIpointer();

  iExit = false;
  iCont = false;
  
  for(int i=0;i<4;i++)
  { 
    brktbl[i] = false; 
  }
  
  noBreakFlag = true;
  verboseFlag = true;
  new G4UIxvtMessenger(this);

  number = 100;
  currentPosition = 0;
  lastPosition = 100;
}


////////////////////////////////////////////////////////////////////////////
//
// ~G4UIxvt - Destructor
//
////////////////////////////////////////////////////////////////////////////
G4UIxvt::~G4UIxvt()
{

}

////////////////////////////////////////////////////////////////////////////
// 
// Notify
//
// Notifies the user of the state changes of the system.
//
////////////////////////////////////////////////////////////////////////////
//G4bool G4UIxvt::Notify(G4ApplicationState requestedState)
//{
//  G4StateManager * statM = G4StateManager::GetStateManager();
//  G4ApplicationState previousState = statM->GetPreviousState();
//
//  if(verboseFlag) 
//  {
//    writeGeantToXvt("========== STATE CHANGE ===========");
//    writeGeantToXvt("FROM: " + statM->GetStateString(previousState));
//    writeGeantToXvt("  TO: " + statM->GetStateString(requestedState));
//    writeGeantToXvt("===================================");
//  }
//
//  if(requestedState==Pause || breakRequested(previousState,requestedState))
//  { 
//    additionalSession(); 
//  }
//  
//  return true;
//}
 
////////////////////////////////////////////////////////////////////////////
// 
// breakRequested
//
// Notifies and puts breaks in the Geant4 run.
//
////////////////////////////////////////////////////////////////////////////
G4bool G4UIxvt::breakRequested(G4ApplicationState previousState,
                                    G4ApplicationState requestedState)
{
  if(noBreakFlag) return false;

  G4bool ans = false;

  if(previousState==Idle && requestedState==GeomClosed && brktbl[0])
  { 
    writeGeantToXvt("====== BREAKPOINT NOTIFICATION ======");
    writeGeantToXvt("BreakPoint at: BEGIN OF RUN");
    writeGeantToXvt("Enter 'continue' command to proceed");
    writeGeantToXvt("=====================================");
    
    ans = true;
  }

  if(previousState==GeomClosed && requestedState==EventProc && brktbl[1])
  { 
    writeGeantToXvt("====== BREAKPOINT NOTIFICATION ======");
    writeGeantToXvt("BreakPoint at: BEGIN OF EVENT");
    writeGeantToXvt("Enter 'continue' command to proceed");
    writeGeantToXvt("=====================================");
    
    ans = true;
  }

  if(previousState==EventProc && requestedState==GeomClosed && brktbl[2])
  {     
    writeGeantToXvt("====== BREAKPOINT NOTIFICATION ======");
    writeGeantToXvt("BreakPoint at: END OF EVENT");
    writeGeantToXvt("Enter 'continue' command to proceed");
    writeGeantToXvt("=====================================");
    
    ans = true;
  }

  if(previousState==GeomClosed && requestedState==Idle && brktbl[3])
  { 
    writeGeantToXvt("====== BREAKPOINT NOTIFICATION ======");
    writeGeantToXvt("BreakPoint at: END OF RUN");    
    writeGeantToXvt("Enter 'continue' command to proceed");
    writeGeantToXvt("=====================================");
    
    ans = true;
  }

  return ans;
}

////////////////////////////////////////////////////////////////////////////
//
// set_verbose
//
// Toggles the Verbose state of the GUI concerning the state machine.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::set_verbose(void)
{
  verboseFlag = !verboseFlag;
  
  if (verboseFlag)   
  {
    writeGeantToXvt("Verbose mode switched ON");
  }
  else
  {
    writeGeantToXvt("Verbose mode switched OFF");
  }
}


////////////////////////////////////////////////////////////////////////////
// 
// SessionStart
//
// Informs user of Geant Kernel starting and does a 'brief' command listing.
// Update: 12/08/97 - now used as main control of Geant (get command is 
//         currently called from here NOT from the UIManager.
//
////////////////////////////////////////////////////////////////////////////
//void G4UIxvt::SessionStart(void)
G4UIsession*  G4UIxvt::SessionStart(void)
{
  openConnections();
  writeGeantToXvt("Connection to Geant open...\0");    
  
  G4UImanager * UI = G4UImanager::GetUIpointer();
  G4UIcommandTree * comTree = UI->GetTree();

  fillArrayEntries(comTree, 0);

  G4cout << " --------- List of Array Entries ---------- " << endl;

  listCommandArray();
  briefListCommands();
  
  // --- Enter main program control --- //
  
  iExit = true;

  while( iExit )
  {
    G4String newCommand = GetCommand();
    if(newCommand.length()>1) UI->ApplyCommand(newCommand);
  }    
  return NULL;
}


////////////////////////////////////////////////////////////////////////////
// 
// additionalSession
//
// 
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::additionalSession(void)
{
  writeGeantToXvt("=== Entering Addtional Session ===");
  iCont = true;
  while( iCont )
  {
    G4String newCommand = GetCommand();
    if(newCommand.length()>1) UI->ApplyCommand(newCommand);
  }
}


////////////////////////////////////////////////////////////////////////////
//
// SessionTerminate
//
// Informs user that Geant has terminated.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::SessionTerminate(void)
{
  writeGeantToXvt("XVT <-> Geant4 terminating...");
}


////////////////////////////////////////////////////////////////////////////
//
// openConnections                 
//
// This method opens the IPC connection.   
// 
// Open the named pipes in the current directory using the 
// appropriate permissions.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::openConnections(void)
{          
  int o_flags;
  
  o_flags = O_WRONLY | O_NONBLOCK;
  
  if((fd_GeantToXvt = open(GeantToXvt, o_flags)) == -1)
  {
    errorHandler("GeantToXvt pipe opening");
  }
  
  o_flags = O_RDONLY | O_NONBLOCK;
  
  if((fd_XvtToGeant = open(XvtToGeant, o_flags)) == -1)
  {
    errorHandler("XvtToGeant pipe opening"); 
  } 
}


////////////////////////////////////////////////////////////////////////////
//                   
// checkXvtToGeantPipe
//
// This method checks to see if the XvtToGeant pipe contains anything.
// Thus meaning it will be readable without blocking.
//
// Returns 0 if nothing present, 1 otherwise. 
// 
////////////////////////////////////////////////////////////////////////////
int G4UIxvt::checkXvtToGeantPipe(void)
{
  int nbytes;
  ioctl(fd_XvtToGeant, FIONREAD, &nbytes);
  
  if (nbytes > 0)
    return 1;
  else
    return 0;
}


////////////////////////////////////////////////////////////////////////////
//                   
// readXvtToGeant
//
// This method grabs the data in the XvtToGeant named pipe, and 
// returns it to the caller.
//
////////////////////////////////////////////////////////////////////////////
G4String G4UIxvt::readXvtToGeant(void)
{
  int bytes_read = read(fd_XvtToGeant, textDump, number);
  textDump[bytes_read] = '\0'; 
  G4String newString = textDump;
  return newString;
}
 
 
////////////////////////////////////////////////////////////////////////////
//
// writeGeantToXvt                 
//
// This method writes the passed string to the GeantToXvt named pipe.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::writeGeantToXvt(const G4String& theString)
{
  int write_check = write(fd_GeantToXvt, theString, number);
    
  if(write_check == -1)
    G4cout << "ERROR WITH WRITE!!" << endl;
}


///////////////////////////////////////////////////////////////////////////////
//
// errorHandler
//
// A general purpose error handler.  It can be called by any method.  A string
// is passed so that this method simply prints a message via G4cout that explains 
// what the error is.
//
///////////////////////////////////////////////////////////////////////////////
void G4UIxvt::errorHandler(const G4String& theError)
{
  G4cout << "Error with " << theError << endl;  
}


////////////////////////////////////////////////////////////////////////////
//
// GetCommand
//
// This method is called by SessionStart to get a command from the user.
//
// It enters a loop to repeatedly check the Geant end of the named pipe.
// When something is present in the pipe (sent by XVT GUI) it then grabs it 
// and returns the string to the caller.
//
////////////////////////////////////////////////////////////////////////////
G4String G4UIxvt::GetCommand()
{
  G4String newCommand('\0');
  G4UImanager * UI = G4UImanager::GetUIpointer();
  while( 1 )
  {
    int Done = 0;
  
    while (!Done)
    {
      if (checkXvtToGeantPipe())
      {
        newCommand = readXvtToGeant();
        Done = 1;
      }
      else
        ; // Pipe is empty so loop again
    }

    G4String nC = newCommand;
    
    if(nC == "getCommands")
    {
      sendArrayEntriesToXvt();
    }
    else if( nC == "/control/exit" )
    { 
      if( iCont )
      { 
        writeGeantToXvt("You are now processing RUN.");
        writeGeantToXvt("Please abort it using \"/run/abort\" command first");
        writeGeantToXvt("and use \"continue\" command until the application");
        writeGeantToXvt("becomes Idle.");        
      }
      else
      {
        iExit = false;
        newCommand = "/";
        break;
      }
    }    
    else if( nC(0,4) == "cont" )
    { 
      iCont = false;
      newCommand = "/";
      break;
    }    
    else
    { 
      break; 
    }
  }

  if( newCommand(0) != '/' )
    G4cerr << "IT GOT THERE!!!!" << endl;
    // If you ever see the error message in the output something has gone
    // wrong.

  return newCommand;

}


////////////////////////////////////////////////////////////////////////////
//
// codeGen
//
// Traverses the command tree and retrieves the directories and 
// subdirectories of available commands.
// 
// This does a 'brief' listing, i.e. does not List the guidance or 
// parameters.
//
// I have updated this method from the TCL implementation for XVT.
// The differences are that G4cout is not used and some string concatenation
// is Done differently (i.e. not with G4cout).
//
// This is still however not really needed and might be deleted.
// 
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::codeGen(G4UIcommandTree * tree, int level)
{ 
  int treeEntry;
  int commandEntry; 
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();
  G4String commandPath;
  G4String pathName; //tree name
  G4UIcommandTree * t;

  // --- this loop lists the current directory and its commands --- //
  
  for(int com=0; com<commandEntry; com++)
  {
    commandPath = tree->GetCommand(com+1)->GetCommandPath();
    if(level==0) 
    {
      // This is not normally executed    
      writeGeantToXvt(commandPath);  // Unsure as to what this does         
    }
    else     
    {      
      writeGeantToXvt(commandPath + " command");  // this gives the commands
    }
  }

  if(treeEntry == 0) 
    return; //end recursion

  for(int i=0; i< treeEntry; i++)
  {
    t = tree->GetTree(i+1);
    pathName =  t->GetPathName();   
    if(level == 0) 
    {
      writeGeantToXvt(pathName);  // This is the directory name
    }
    else
    {      
      writeGeantToXvt(pathName + " cascade"); // This is the sub dirs      
    }
      
    codeGen(t, level+1);
  }
}


////////////////////////////////////////////////////////////////////////////
//
// briefListCommands
//
// This method uses codeGen to traverse the tree and List the commands
// and directories - it does a brief listing, i.e. no parameters or guidance
// accompanying the text.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::briefListCommands(void)
{
  writeGeantToXvt("--- Brief List of available commands ---");
  G4UIcommandTree * comTree = UI->GetTree();
  codeGen(comTree, 0);
}


////////////////////////////////////////////////////////////////////////////
//
// ListCurrentDirectory
//
// This is the XVT specific implementation of the 'List Current' method
// from the G4UICommandTree class.
// 
// It is almost the same only instead of any couts all the data is written
// to the named pipe for XVT to process.
//
// Although this is entirely XVT specific it may be deleted because it is
// not really needed. (most directory navigation methods that were so
// important in the terminal session aren't needed in this XVT GUI
// implementation.)
//
// I have left it here in case similar code is needed in the future - this
// demonstrates how to access the command tree etc.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::listCurrentDirectory(void)
{
  G4UIcommandTree * comTree = UI->GetTree();
  G4String path = comTree->GetPathName();
  
  writeGeantToXvt("Command Directory Path: " + path);

  // --- Guidance code --- //

  G4UIcommand * tempGuidance = (G4UIcommand *)comTree->GetGuidance();
 
  if (tempGuidance != NULL)
  {
    writeGeantToXvt("Command: " + tempGuidance->GetCommandPath());
    writeGeantToXvt("Guidance: ");

    int n_guidanceEntry = tempGuidance->GetGuidanceEntries();

    for (int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ )
    {
      writeGeantToXvt(tempGuidance->GetGuidanceLine(i_thGuidance));
    }

    int n_parameterEntry = tempGuidance->GetParameterEntries();
  
    if(n_parameterEntry > 0)
    {
      for( int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ )
      {
        G4UIparameter * tempParameter = (G4UIparameter *)tempGuidance->GetParameter(i_thParameter);
    
        writeGeantToXvt("Parameter: " + tempParameter->GetParameterName());
        writeGeantToXvt(tempParameter->GetParameterGuidance());
        writeGeantToXvt("Parameter type: " + tempParameter->GetParameterType());
      
        // --- remember to fix the below bug --- //
      
        // writeGeantToXvt("Omittable: " + tempParameter->IsOmittable());
      
        writeGeantToXvt("Default Value: " + tempParameter->GetDefaultValue());
        writeGeantToXvt("Parameter Range: " + tempParameter->GetParameterRange());
      }
    }                    
  }
  
  // --- End of Guidance code --- //
  
  writeGeantToXvt("Sub Directories: ");
  
  int n_treeEntry = comTree->GetTreeEntry();
  
  for( int i_thTree = 0; i_thTree < n_treeEntry; i_thTree++ )
  {
    G4UIcommandTree * tempTree = comTree->GetTree(i_thTree + 1);    
    writeGeantToXvt(tempTree->GetPathName() + tempTree->GetTitle());
  }                     
  
  writeGeantToXvt("Commands : ");
  
  int n_commandEntry = comTree->GetCommandEntry();
  
  for( int i_thCommand = 0; i_thCommand < n_commandEntry; i_thCommand++ )
  {
    G4UIcommand * tempCommand = comTree->GetCommand(i_thCommand + 1);    
    writeGeantToXvt(tempCommand->GetCommandName() + tempCommand->GetTitle());
  }                     
}

////////////////////////////////////////////////////////////////////////////
// Methods for handling the command information                           //
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//
// sendArrayEntriesToXvt
//
// This goes through the commandArray when called and sends all the command
// information down the named pipe to XVT.
// 
// It must first send the control command 'comTree' to alert XVT of the 
// following commands - XVT can then take action to fill its own 'commandArray'
// with the ensuing strings passed to it.
//
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::sendArrayEntriesToXvt(void)
{
  writeGeantToXvt("comTree");
  
  char theInt[5];
  sprintf(theInt, "%d", currentPosition);
  
  writeGeantToXvt(theInt);

  // --- Output all the Data from the array to the named pipe --- //
  
  for (int i = 0; i < currentPosition; i++)
  {
    // --- Send flag to XVT --- //
    sprintf(theInt, "%d",commandArray[i].flag);     
    writeGeantToXvt(theInt);
    
    // --- Send name and guidance --- //
    writeGeantToXvt(commandArray[i].name);
    writeGeantToXvt(commandArray[i].guidance);
    
    // --- Send number of parameters --- //
    sprintf(theInt, "%d", commandArray[i].numOfParameters);
    writeGeantToXvt(theInt);
    
    // --- if there are parameters send them next --- //
    if (commandArray[i].numOfParameters > 0)
    {
      for(int x = 0; x < commandArray[i].numOfParameters; x++)
      {
        writeGeantToXvt(commandArray[i].parameters[x].name);
        writeGeantToXvt(commandArray[i].parameters[x].guidance);
        writeGeantToXvt(commandArray[i].parameters[x].defaultValue);
        writeGeantToXvt(commandArray[i].parameters[x].range);
        writeGeantToXvt(commandArray[i].parameters[x].candidate);
        writeGeantToXvt(commandArray[i].parameters[x].type);
        writeGeantToXvt(commandArray[i].parameters[x].omittable);      
      }    
    }             
    
    // --- Send number --- //
    sprintf(theInt, "%d", commandArray[i].commandNumber);    
    writeGeantToXvt(theInt);
  }   
}


////////////////////////////////////////////////////////////////////////////
//
// fillArrayEntries
//
// This traverses the command tree and fills in the command structures in 
// the class variable 'commandArray'.
//
// Note: for flag value:
//
//       0 = Directory
//       1 = command
//       2 = subdirectory
// 
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::fillArrayEntries(G4UIcommandTree * tree, int level)
{
  int treeEntry;
  int commandEntry; 

  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  G4String commandPath;
  G4String pathName; //tree name
  G4UIcommandTree * t;

  // New variables
  G4UIcommand * tempCommand;
  G4UIparameter * tempParameter;
  
  // --- this loop lists the current directory and its commands --- //
  
  for(int com=0; com<commandEntry; com++)
  {
    tempCommand = (G4UIcommand *)tree->GetCommand(com+1);    
    
    commandPath = tree->GetCommand(com+1)->GetCommandPath();        
    if(level==0) 
    {
      // --- This code shouldn't ever get executed --- //
      G4cout << "Printing command Path" << commandPath << endl;
    }
    else     
    {
      commandArray[currentPosition].flag = 1;
      commandArray[currentPosition].name = commandPath;      
      commandArray[currentPosition].guidance = tempCommand->GetTitle();      
      commandArray[currentPosition].numOfParameters = tempCommand->GetParameterEntries();
            
      int numParam = tempCommand->GetParameterEntries();
      
      if(numParam != 0)
      {              
        for(int x = 0; x < numParam; x++)
        {
          tempParameter = (G4UIparameter*)tempCommand->GetParameter(x); 

          commandArray[currentPosition].parameters[x].name = tempParameter->GetParameterName();
          commandArray[currentPosition].parameters[x].guidance = tempParameter->GetParameterGuidance();
          commandArray[currentPosition].parameters[x].defaultValue = tempParameter->GetDefaultValue();
          commandArray[currentPosition].parameters[x].range = tempParameter->GetParameterRange();
          commandArray[currentPosition].parameters[x].candidate = tempParameter->GetParameterCandidate();
          commandArray[currentPosition].parameters[x].type = tempParameter->GetParameterType();
//          commandArray[currentPosition].parameters[x].omittable = tempParameter->IsOmittable();
        }      
      }
      
      commandArray[currentPosition].commandNumber = currentPosition;

      currentPosition = currentPosition + 1;
    }
  }

  if(treeEntry == 0) 
    return; //end recursion
  
  G4UIcommand * treeGuidance;  
  
  for(int i=0; i< treeEntry; i++)
  {
    t = tree->GetTree(i+1);
    pathName =  t->GetPathName();   
    
    treeGuidance = (G4UIcommand *)t->GetGuidance();
    
    if(level == 0) 
    {
      commandArray[currentPosition].flag = 0;
      commandArray[currentPosition].name = pathName;
      commandArray[currentPosition].guidance = treeGuidance->GetTitle();      
      commandArray[currentPosition].numOfParameters = 0;
      commandArray[currentPosition].commandNumber = currentPosition;

      currentPosition = currentPosition + 1;
    }
    else
    {      
      commandArray[currentPosition].flag = 2;
      commandArray[currentPosition].name = pathName;
      commandArray[currentPosition].guidance = treeGuidance->GetTitle();      
      commandArray[currentPosition].numOfParameters = 0;
      commandArray[currentPosition].commandNumber = currentPosition;
      
      currentPosition = currentPosition + 1;
    }
      
    fillArrayEntries(t, level+1);
  }
}


////////////////////////////////////////////////////////////////////////////
// A debugging routine                                                    //
//                                                                        //
// This method is purely used for debugging purposes - will be deleted    //
// from the final code (when we are sure it is bug free)                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
void G4UIxvt::listCommandArray(void)
{
  for(int i = 0; i < currentPosition; i++)
  {
    G4cout << "CommandArray, Entry-> " << i << endl
         << "Flag-> " << commandArray[i].flag << endl
         << "Name-> " << commandArray[i].name << endl
         << "Guidance-> " << commandArray[i].guidance << endl
         << "Number of Parameters-> " << commandArray[i].numOfParameters << endl;
         
         int numParam = commandArray[i].numOfParameters;
         
         if (numParam != 0)
         {
           for(int x = 0; x < numParam; x++)
           {
             G4cout << "Parameters for this command: " << endl << endl
                  << "Parameter Name-> "          << commandArray[i].parameters[x].name         << endl
                  << "Parameter Guidance-> "      << commandArray[i].parameters[x].guidance     << endl
                  << "Parameter Default Value-> " << commandArray[i].parameters[x].defaultValue << endl
                  << "Parameter Range-> "         << commandArray[i].parameters[x].range        << endl
                  << "Parameter Candidate-> "     << commandArray[i].parameters[x].candidate    << endl
                  << "Parameter Type-> "          << commandArray[i].parameters[x].type         << endl
                  << "Parameter Omittable-> "     << commandArray[i].parameters[x].omittable    << endl 
                  << endl << endl;
           }
         }
         
         G4cout << "Number-> " << commandArray[i].commandNumber << endl << endl;
  }
}

#endif
