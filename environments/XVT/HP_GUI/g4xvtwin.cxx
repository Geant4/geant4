// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: g4xvtwin.cxx,v 1.1 1999-01-07 16:05:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/****************************************************************************
 *
    Copyright (c) 1996

    File Name:    G4XvtWin.cxx

    Implements a window class to display data for a derived document class.

    The properties of this window were specified using XVT-Architect, and the
    window is created at runtime by requesting an instance from the factory.
 *
 ***************************************************************************/

#include "XVTPwr.h"
#include "AppDef.h"

#include G4XvtWin_i
#include CNavigator_i
//#include CMenuBar_i

#include NScrollText_i
#include CStatusBar_i

#include NLineText_i
#include NEditControl_i

// Define the Run-Time Type Identification for this class
// PWRRegisterClass1( G4XvtWin, G4XvtWinID, "G4XvtWin", CWindow )

#include CToolBarAttachment_i
#include CAttachment_i
#include CAttachmentFrame_i
#include G4XvtGUI_i
#include CToolPalette_i
#include CViewSink_i
#include CTaskWin_i
#include CMouseManager_i
#include CImage_i
#include CScroller_i
#include CMenuButton_i

#include "faccmd.h"

///////////////////////////////////////////////////////////////////////////////
//
// Factory Window class constructor
//
// Instantiate a window object with properties known to the factory
//
///////////////////////////////////////////////////////////////////////////////
G4XvtWin::G4XvtWin(
	CDocument *      theDocument,
	const CRect&     theRegion,
	const CStringRW& theTitle,
	long             theAttributes,
	WIN_TYPE         theWindowType,
	int              theMenuBar,
	WINDOW           theParentWindow )
	: CWindow( theDocument, theRegion, theTitle, theAttributes,
		theWindowType, theMenuBar, theParentWindow )
{
	// Ask the factory to create the contents of this window.
	// Create all view objects that had the AutoCreate option
	// specified in XVT-Architect.
	SuspendSizing();
	GUIFactory.DoCreateViews( this, G4Win, FALSE, &itsData );

	// Add code here to Initialize the views that are
	// nested in this window via the factory definition.

	// You may want to put initialization code in the I4XvtWin()
	// method, and then explicitly call the method any time you
	// want the window to be set to its initial state               
        
        ///////////////////////////////////////////////////////////////////////
        // Carry out initialisation - order IS important                     //
        ///////////////////////////////////////////////////////////////////////
        
        initIndex();
        initStatusBar();
        initTimer();
        initTextWindows();
        initIPC();
        initFrame();              
        initCommandPalette();
        
        // Arrange all attachments in the frame so they fit nicely
        itsAttachmentFrame->ArrangeAttachments();

	// Ensure objects created in the constructor are updated.
	ResumeSizing();
	GetNavigator()->InitFocus();
	xvt_dwin_invalidate_rect( GetXVTWindow( ), NULL );
}


///////////////////////////////////////////////////////////////////////////////
//
// Window class destructor
//
// Delete the internal data objects for this class
//
///////////////////////////////////////////////////////////////////////////////
G4XvtWin::~G4XvtWin( void )
{
  G4cout << "Cleaning up...." << endl;

  xvt_timer_destroy(theTimer);  
  destroyNamedPipes();  
  destroyCommandPalette();  
  destroyDragSource();
}


///////////////////////////////////////////////////////////////////////////////
//
// Factory Window class initializer
//
// Sets the properties of the window and any nested views to a standard state
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtWin::I4XvtWin( void )
{
	// Add any code to Initialize the window here.  This method is
	// here for convenience only, and it is not called automatically.
	// Your code must explicitly call this method if you put code here.

	return( TRUE );
}


///////////////////////////////////////////////////////////////////////////////
//
// DoCommand
//
// Handle commands here that have to do with the display
// of the data, or that create, destroy or otherwise manage the
// nested views that the window uses to help display the data.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::DoCommand( long theCommand, void * theData )
{
	switch( theCommand )
	{
	// --- Button descriptions --- //
	
	case START_INcmd:
	  setStatus("Start Geant4...");
	  break;
	
	case START_OUTcmd:
	  resetStatus();
	  break;
	  
	case STOP_INcmd:
	  setStatus("Stop Geant4...");
	  break;
	
	case STOP_OUTcmd:
	  resetStatus();
	  break;
	
	case PAUSE_INcmd:
	  setStatus("Pause Geant4...");
	  break;
	
	case PAUSE_OUTcmd:
	  resetStatus();
	  break;
	
	case RESUME_INcmd:
	  setStatus("Resume from Pause...");
	  break;
	
	case RESUME_OUTcmd:
	  resetStatus();
	  break;
	
	case QUIT_INcmd:
	  setStatus("Quit Geant4...");
	  break;
	
	case QUIT_OUTcmd:
	  resetStatus();
	  break; 

        case SAVE_LOG_INcmd:
	  setStatus("Save Log to file...");
	  break;
	
	case SAVE_LOG_OUTcmd:
	  resetStatus();
	  break;
	
	case SAVE_AS_INcmd:
	  setStatus("Save Log AS...");
	  break;
	  
	case SAVE_AS_OUTcmd:
	  resetStatus();
	  break;
	
	case PRINT_INcmd:
	  setStatus("Print Log...");
	  break;
	  
	case PRINT_OUTcmd:
	  resetStatus();
	  break;
	  
	case HELP_INcmd:
	  setStatus("Get Help on currently selected command...");
	  break;
	  
	case HELP_OUTcmd:
	  resetStatus();
	  break;

        case EXECUTE_INcmd:
          setStatus("Execute currently selected command...");
          break;
          
        case EXECUTE_OUTcmd:
          resetStatus();
          break;
        
        case GET_COMMANDS_INcmd:
          setStatus("Retrieve all the Geant 4 available commands...");
          break;
          
        case GET_COMMANDS_OUTcmd:
          resetStatus();
          break;          
        
        case COM_INcmd:
          setStatus("Enter command through text box...");
          break;
          
        case COM_OUTcmd:
          resetStatus();
          break;          
        
	case NULLcmd:
		// stop propagation of "undefined" commands
	case NATIVEViewCmd:
		// stop commands only of interest to CScroller
		break;

	default:        
		CWindow::DoCommand( theCommand, theData );
		break;
	
	}
}


///////////////////////////////////////////////////////////////////////////////
//
// DoMenuCommand
//
// Handle menu items here that have to do with the display
// of the data, or that create, destroy or otherwise manage the
// nested views that the window uses to help display the data.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::DoMenuCommand(
	MENU_TAG theMenuItem,
	BOOLEAN  isShiftKey,
	BOOLEAN  isControlKey )
{
	switch( theMenuItem )
	{
	case M_Retrieve_Commands:
	  G4cout << "Retrieving commands" << endl;
	  writeData("getCommands");
	  G4cout << "Commands retrieved" << endl;	
	  
	  // --- Now prevent user from clicking that button again --- //
	  // --- by disabling it and redrawing it greyed out      --- //
	  
	  itsData.itsgetComButton->Disable();
	  itsData.itsgetComButton->DoDraw(itsData.itsgetComButton->GetGlobalFrame());

	  break;
	  		  
	case M_Commands:
	  regenerateCommandPalette();
	  G4cout << "Would Retrieve Commands Palette if this was selected"
	       << endl;
	  break;	

	case STARTcmd:
	case M_Start_Geant:
          if(!G4Running)
          {
            G4cout << "Calling startSlave" << endl;
            destroyNamedPipes();
            createNamedPipes();
            startSlave();
          }
	  break;
	  
	case M_Stop_Geant:
	case STOPcmd:
	  break;
	
	case M_Pause_Geant:
	case PAUSEcmd:
	  break;	
	
	case M_Resume_Geant:
	case RESUMEcmd:
	  break;
        
        case M_Quit_Geant:
        case QUITcmd:
          G4cout << "Quit button pressed" << endl;
          writeData("/control/exit");
	  G4Running = 0;
          break;
          
	case M_On_command:
	case HELPcmd:
	  showHelp();	  
	  break;
        
        case EXECUTEcmd:
        {
          G4cout << "Command would be executed when this was pressed" << endl;
          ExecuteCommand();                              
        }
        break;

        case COMcmd:
        {
          CStringRW tempString = getTextCommand();
          if(tempString != "nothing")
            writeData(tempString);          
        }  
        break;
        
	case NULLcmd:
		// stop propagation of "undefined" commands
		break;

	default:
		CWindow::DoMenuCommand( theMenuItem, isShiftKey, isControlKey );
		break;
	}
}


///////////////////////////////////////////////////////////////////////////////
//
// UpdateMenus
//
// Respond to notifications that the window just gained input focus
// or that the "dirty bit" of the data class just changed state.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::UpdateMenus( CMenuBar * /* theMenuBar */ )
{
	// Two standard behaviors to implement:
	//   1) Enable the File|Save menu option only when the data needs
	//      to be saved
	//   2) Enable the Edit|Paste menu option if this window supports
	//      pasting and there is something paste-able on the clipboard

	// theMenuBar->SetEnabled( M_FILE_SAVE, itsDocument -> NeedsSaving( ) );
	// theMenuBar->SetEnabled( M_EDIT_PASTE, ... );
}


///////////////////////////////////////////////////////////////////////////////
//
// DoUpdateModel
//
// Respond to notifications that the model has changed for any
// data this window is registered as a dependent of.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::DoUpdateModel(
	long           /* theControllerId */,
	long           theCommand,
	const CModel * /* theModel */ )
{
	switch( theCommand )
	{
	case NULLcmd:    // stop propagation of "undefined" commands
		break;

	default:
		xvt_dm_post_note( "Command (%ld) not handled in dependent window %s",
			theCommand, GetTitle( ).data( ) );
		break;
	}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////// Initialisation methods (called from the constructor on startup) ////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initIndex                                                                 //
//                                                                           //
// Initialise the command base indexing variables                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initIndex(void)
{
  currentPosition = 0;
  lastPosition = 100;
  reading = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initStatusBar                                                             //
//                                                                           //
// Initialise the status bar                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initStatusBar(void)
{
  CStatusBar* SBar = itsData.itsprogStatus;

  SBar->AppendField(NULL_FIELD,"  ", 30);
  SBar->AppendField(STATUS_FIELD, "Status: initialisation");
  SBar->AppendField(COORDS_FIELD, "0, 0", 100);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initTimer                                                                 //
//                                                                           //
// Set up the timer.                                                         //
//                                                                           //
// Alter interval value below to change the time interval between calls to   //
// the DoTimer method                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initTimer(void)
{
  interval = 1000; 
  theWin = (itsData.itsmainTextArea->GetCWindow())->GetXVTWindow();
  theTimer = xvt_timer_create(theWin, interval);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initIPC                                                                   //
//                                                                           //
// Set the value of the pipe reading variable (i.e. how many bytes are read  //
// in a single read operation and then create the named pipes and start      //
// Geant.                                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initIPC(void)
{
  number = 100;
  createNamedPipes();
  startSlave();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initFrame                                                                 //
//                                                                           //
// Create an attachment frame and allow attachment to all sides.             // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initFrame(void)
{
  itsAttachmentFrame = new CAttachmentFrame( itsData.itsDropArea,
                                             itsData.itsDropArea->GetLocalFrame(), 
                                             CAttachmentFrame::ATTACH_ALL );
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////// My methods for IPC etc //////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// createNamedPipes
//
// This method creates the named pipes required in the application.
//  
// It first creates the pipes and then opens them.  
//
// The XvtToGeant pipe with write permission only, 
// the GeantToXvt with read permission only.
//
// The program uses the errorHandler method to report errors because the
// program should Terminate if IPC setup fails.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::createNamedPipes(void)
{
  int o_flags;
  
  if(mkfifo(XvtToGeant, FILE_MODE) == -1)
  {
    errorHandler("XvtToGeant pipe creation");
  }    
  
  o_flags = O_WRONLY | O_NONBLOCK;
  
  if((fd_XvtToGeant = open(XvtToGeant, o_flags)) == -1)
  {
    errorHandler("XvtToGeant pipe opening");
  }
  
  if(mkfifo(GeantToXvt, FILE_MODE) == -1)
  {
    errorHandler("GeantToXvt pipe creation");
  }

  o_flags = O_RDONLY | O_NONBLOCK;
  
  if((fd_GeantToXvt = open(GeantToXvt, o_flags)) == -1)
  {
    errorHandler("GeantToXvt pipe opening"); 
  }

}


///////////////////////////////////////////////////////////////////////////////
//
// destroyNamedPipes
//
// This method closes and unlinks all the open named pipes.
//  
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::destroyNamedPipes(void)
{
  close(fd_XvtToGeant);
  close(fd_GeantToXvt);
  
  unlink(XvtToGeant);
  unlink(GeantToXvt);
}

      
///////////////////////////////////////////////////////////////////////////////
//
// startSlave
//
// This method forks and exec's the other process.
// 
// This method also uses the errorHandler to report errors because if the 
// fork or Exec fail the program should Terminate.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::startSlave(void)
{
  pid_t pid;
  
  if((pid = fork()) < 0)
    errorHandler("Fork");
    
  else if(pid > 0)
       {
         G4Running = 1; // --- set flag --- //
       }
       
       else      // in child process
       {         
         if(execlp("xvt.benchmark", "", (char *) 0) < 0)
           errorHandler("Exec");
       }
}


///////////////////////////////////////////////////////////////////////////////
//
// DoTimer
//
// Every time the timer is activated this method
// checks to see if the file descriptor can be read without blocking (i.e.
// there is something in the pipe), if there is it grabs the data from the 
// pipe, if you are not in the process of reading all the commands from Geant 
// it then checks if the string is equal to "comTree" which indicates all the
// commands are coming from Geant and different action should be taken.  If it
// isn't comTree then it simply displays the string in the Main log Area.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::DoTimer(long theTimerId)
{
  int nbytes;
  ioctl(fd_GeantToXvt, FIONREAD, &nbytes);
  
  if (nbytes > 0 && !reading)
  {
    int bytes_read = read(fd_GeantToXvt, textDump, number);
    textDump[bytes_read] = '\0';
    
    CStringRW tempString = textDump;
    
    if (tempString == "comTree")
    {
      reading = 1; 
      fillLocalArray();
      buildCommandPalettes();
      reading = 0;
      listCommandArray();	  
      activateButtons();
    }
    else
    {
      itsData.itsmainTextArea->Append(textDump);
    }  
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// errorHandler
//
// A general purpose error handler.  It can be called by any method.  A string
// is passed so that this method can post a dialog box to the user stating 
// where the program fault occurs.
//
// KEY POINT: This method terminates execution of the program when called.  If
//            the error that needs reporting is not a 'critical' error - use
//            G4cout.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::errorHandler(CStringRW theError)
{
  xvt_dm_post_note("Error with %s", theError);
  xvt_app_destroy();
}


///////////////////////////////////////////////////////////////////////////////
//
// writeData
//
// Writes the passed string to the Command log window and the XvtToGeant 
// pipe - i.e. will send the string to Geant (commands!).
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::writeData(CStringRW theString)
{
  itsData.itsCommandLog->Append(theString);

  int write_check = write(fd_XvtToGeant, theString, number);
    
  if(write_check == -1)
    G4cout << "ERROR WITH WRITE!!" << endl;
}

///////////////////////////////////////////////////////////////////////////////
//
// initCommandPalette
//
// Basically just allocate storage space for a tool palette.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initCommandPalette(void)
{
  commandPalette = new CToolPalette("Commands"); 
}

///////////////////////////////////////////////////////////////////////////////
//
// destroyCommandPalette
//
// Free up storage space allocated for the command palette.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::destroyCommandPalette(void)
{
  delete(commandPalette); 
}

///////////////////////////////////////////////////////////////////////////////
//
// buildCommandPalettes
//
// This method builds the palettes of commands available to the user. 
// 
// This basically parses the commandArray that was received from Geant and
// builds up the command palette 'on-the-fly' from this information.  
//
// It checks the flag of the command in the array to see whether it is
// a 'directory', a command or a cascade.
//
// Once it has this information it can decide how to handle it.. If it is a
// 'directory' - i.e. a category of command it must have its own sub-palette
// so this can be constructed.  If it is a command it must be a button in the
// current sub-palette and if it is a cascade it must be a sub-palette of the
// current palette - see code for details.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::buildCommandPalettes(void)
{
  CToolPalette* tempPalette;
  CToolPalette* currentParent;

  for(int i = 0; i < currentPosition; i++)
  {
    switch(commandArray[i].flag)
    {
      case 0:
        tempPalette = new CToolPalette( commandArray[i].name );
        commandPalette->AppendTool( new CImage("temp.bmp"), commandArray[i].commandNumber, tempPalette);
        currentParent = tempPalette;
        break;

      case 1:
        currentParent->AppendTool(new CImage("red.bmp"), commandArray[i].commandNumber, TRUE);
        break;
      
      case 2:
        tempPalette = new CToolPalette( commandArray[i].name );
        currentParent->AppendTool( new CImage("oval.bmp"), commandArray[i].commandNumber, tempPalette);        
        currentParent = tempPalette;
        break;
    }
  }                           
  
  // --- set up drag and drop etc and where the palette is initially --- //
  // --- attached.                                                   --- //

  // --- Create a visual representation of the palette. The Popup method --- //
  // --- creates a popup window for the palette and returns a pointer to --- //
  // --- the attachment which will manage the palette.                   --- //
  
  CAttachment* commandAttachment = commandPalette->Popup( 
                this, 
                CPoint(), 
                NULLcmd, 
                TRUE, 
                itsAttachmentFrame );

  // --- Attach the attachment to the left side of the window --- //
  // --- I have left this line commented out as I like the palette to 
  // --- appear popped up instead of attached.
  
  // commandAttachment->Attach( itsAttachmentFrame, CAttachmentFrame::ATTACH_LEFT);

  // ------ Mouse Handler ------- //
  
  commandSource = new CDragSource;
  commandPalette->SetDragSource( commandSource );
  GetMouseManager()->RegisterMouseHandler( commandSource );
  
  // --- Safe Guard against user trying to execute an unselected command by 
  // --- having a pre-set one - this stops core dumps when the program tries
  // --- to detect what tool is selected if one hasn't been
  
  commandPalette->SelectTool(21);  
} 

///////////////////////////////////////////////////////////////////////////////
//
// destroyDragSource
//
// Destroy the space allocated for the drag source. 
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::destroyDragSource(void)
{
  delete(commandSource);
}

  
// --- Status bar code --- //

///////////////////////////////////////////////////////////////////////////////
//
// setStatus
//
// This method places the passed string into the status field of the status 
// bar.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::setStatus(CStringRW theText)
{
  // --- update the status bars text --- //  
  itsData.itsprogStatus->SetFieldStatus(STATUS_FIELD, theText);
}


///////////////////////////////////////////////////////////////////////////////
//
// resetStatus
//
// This method completely blanks the status field of the status bar.
// 
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::resetStatus(void)
{
  // --- reset the fields status to blank --- //
  itsData.itsprogStatus->SetFieldStatus(STATUS_FIELD, " ");
}


// --- Command retrieval methods --- //

///////////////////////////////////////////////////////////////////////////////
//
// retrieveCommand
//
// This method checks the command palette to see which button has been pressed 
// then it will scan the command array to find the command that corresponds to
// the users selection and returns the command string to be sent to Geant.
// 
///////////////////////////////////////////////////////////////////////////////
CStringRW G4XvtWin::retrieveCommand(void)
{
  int selectedTool = commandPalette->GetSelectedTool();
  
  if(selectedTool == NULL)
  {
    xvt_dm_post_note("Command not selected - Nothing executed");  
    G4cout << "Nothing selected" << endl;
  }
  else
  {
    if(commandArray[selectedTool].numOfParameters > 0)
    {
      G4cout << "Inside the if  (num of params > 0)" << endl;
      CStringRW theParameters = getParameters();
      return((commandArray[selectedTool].name) + (" " + theParameters));
    }
    else
    { 
      return commandArray[selectedTool].name;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// initTextWindows
//
// This method simply puts a title in the Text windows to start things off.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::initTextWindows(void)
{
  itsData.itsmainTextArea->Append("-------- Geant 4 Output Log --------");
  itsData.itsCommandLog->Append("-------- Command Log --------");
}

///////////////////////////////////////////////////////////////////////////////
//
// fillLocalArray 
//
// *** Dangerous method - no checks on data, assumes all information is sent
// *** in the correct order from Geant.  
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::fillLocalArray(void)
{
    CStringRW tempString;

    // --- First thing sent is an integer which gives top array position --- //

    tempString = getString();    
    sscanf(tempString, "%d", &currentPosition);
    
    for (int i = 0; i < currentPosition; i++)
    {
      int tempInt;
      
      // --- Grab the flag data member --- //

      tempString = getString();
      sscanf(tempString, "%d", &tempInt);      
      commandArray[i].flag = tempInt;
      

      // --- Grab the name data member --- //

      tempString = getString();      
      commandArray[i].name = tempString;
      

      // --- Grab the guidance data member --- //

      tempString = getString();      
      commandArray[i].guidance = tempString;
      
      // --- Grab the number of parameters data member --- //
      
      tempString = getString();      
      sscanf(tempString, "%d", &tempInt);      
      commandArray[i].numOfParameters = tempInt;     
      
      // --- Grab parameters if they are present --- //
      
      if(tempInt > 0)
      {              
        for(int x = 0; x < tempInt; x++)
        {
          tempString = getString(); 
          commandArray[i].parameters[x].name = tempString;

          tempString = getString();
          commandArray[i].parameters[x].guidance = tempString;

          tempString = getString();
          commandArray[i].parameters[x].defaultValue = tempString;

          tempString = getString();
          commandArray[i].parameters[x].range = tempString;

          tempString = getString();
          commandArray[i].parameters[x].candidate = tempString;

          tempString = getString();
          commandArray[i].parameters[x].type = tempString[0];

          tempString = getString();
          commandArray[i].parameters[x].omittable = tempString;
        }      
      }
      
      // --- Grab the commandNumber data member --- //
       
      tempString = getString();
      sscanf(tempString, "%d", &tempInt);      
      commandArray[i].commandNumber = tempInt;            
    }
}

// --- A debugging routine --- //
void G4XvtWin::listCommandArray(void)
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
                  << endl;
           }
         }                           
         G4cout << "Number-> " << commandArray[i].commandNumber << endl << endl;
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// regenerateCommandPalette 
//
// This method 'nullifies' the original palette pointer and creates a 'new'
// palette.  This is used when the user has closed a palette and wants to 
// open it again.  
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::regenerateCommandPalette(void)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// showHelp
//
// This method is called when the user clicks on the 'help' button.
// 
// It looks which is the current selected tool and retrieves the guidance for
// that tool - then puts this in the text widget and a pop up window.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::showHelp(void)
{
  int selectedTool = commandPalette->GetSelectedTool();
  
  itsData.itsmainTextArea->Append("---------- HELP: ----------");
  itsData.itsmainTextArea->Append(commandArray[selectedTool].guidance);
  itsData.itsmainTextArea->Append("---------------------------");
}

///////////////////////////////////////////////////////////////////////////////
//
// getParameters
//
// This method retrieves parameters from the user via dialog boxes.
//
// ***** Will return a string eventually *******
// see commented out examples below.
// 
///////////////////////////////////////////////////////////////////////////////
CStringRW G4XvtWin::getParameters(void)
{
  char s[100];
  
  // --- Initialise the array --- //

  for (int i=0; i < 100; i++)
  { 
    s[i] = '\0';
  }
  
  if (xvt_dm_post_string_prompt(
        "Type Parameter(s):", s, sizeof(s)) == NULL)
  {
        xvt_dm_post_warning("You cancelled.");        
        return " ";
  }
        
  else
        return s;

}

///////////////////////////////////////////////////////////////////////////////
//
// getTextCommand
//
// This method retrieves a command in text form from the user via dialog boxes.
//
// ***** Will return a string eventually *******
// see commented out examples below.
// 
///////////////////////////////////////////////////////////////////////////////
CStringRW G4XvtWin::getTextCommand(void)
{
  char s[100];
  
  // --- Initialise the array --- //

  for (int i=0; i < 100; i++)
  { 
    s[i] = '\0';
  }
  
  if (xvt_dm_post_string_prompt(
        "Type Command:", s, sizeof(s)) == NULL)
  {
        xvt_dm_post_warning("You cancelled.");
        return "nothing";
  }      

  else
        return s;

}


///////////////////////////////////////////////////////////////////////////////
//
// ExecuteCommand
//
// When the user clicks on the execute button this method is called.
//
// It retrieves the command (with parameters if necessary) and writes it
// down the named pipe to Geant to deal with.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::ExecuteCommand(void)
{
  CStringRW tempString = retrieveCommand();
  G4cout << "Command Retrieved -> " << tempString << endl;
  
  if (tempString == "/control/exit")
  {
    G4cout << "Exit command detected" << endl;
    G4Running = 0;
  }

  G4cout << "Executing command..." << endl;
  
  writeData(tempString);
}


///////////////////////////////////////////////////////////////////////////////
//
// getString
//
// This method grabs the string that is located in the pipe and returns it 
// to the caller - provided there is something present in the pipe.
// 
///////////////////////////////////////////////////////////////////////////////
CStringRW G4XvtWin::getString(void)
{  
  int nbytes;

  int Done = 0;
  
  while (!Done)
  {
    ioctl(fd_GeantToXvt, FIONREAD, &nbytes);
  
    if (nbytes > 0)
    {
      int bytes_read = read(fd_GeantToXvt, textDump, number);
      textDump[bytes_read] = '\0';
      Done = 1;
    }
  }
  
  CStringRW tempString = textDump;
  
  return tempString;
}


///////////////////////////////////////////////////////////////////////////////
//
// activateButtons
//
// This method enables the 'exec' buttons on the toolbar.  It is called just
// after the user has retrieved the commands and built the widget. 
// 
// It first enables the buttons and then refreshes the screen by redrawing 
// them.
// 
///////////////////////////////////////////////////////////////////////////////
void G4XvtWin::activateButtons(void)
{
  itsData.itsexecButton->Enable();
  itsData.itstextExecButton->Enable();
  itsData.itshelpButton->Enable();
  
  itsData.itsexecButton->DoDraw(itsData.itsexecButton->GetGlobalFrame());
  itsData.itstextExecButton->DoDraw(itsData.itstextExecButton->GetGlobalFrame());
  itsData.itshelpButton->DoDraw(itsData.itshelpButton->GetGlobalFrame());  
}
