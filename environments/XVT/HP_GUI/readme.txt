Written by: Simon Prior
Date:       25/08/97

GUI
----

The purpose of this program is to be the first XVT/Geant4 GUI prototype. 

Basic idea: Construct a prototype GUI which interacts with Geant 4.

Future use: After initial investigation into how XVT can interact etc with
            Geant, seeing what features I can use/exploit, it should be 
            possible for future XVT GUI developers to build on the work 
            done adding functionality and design ideas thus constructing one 
            possible GUI for use with Geant 4.  It should be different in
            design/features offered from other implementations so as to
            promote the use of XVT with Geant4 which will cost money (as XVT is
            not free) rather than the use of say TCL which is free.

Discoveries and implementation specifics : 
            
            For this part I shall go through the files which I have edited and
            explain each part - the idea as well as implementation details.
            
            I will give an abstract view to show my overall idea and then a
            detailed one which will highlight implementation specifics.
            
            The two files that do all the work regarding connection to Geant
            and handling of all user input are "g4xvtwin.h" and "g4xvtwin.cxx"
            - the 'window' files.  These are both well commented.                        
            
            Some of the other files have been added to (mainly for the floating
            dockable palettes - see my Geant4 web page or palette documentation
            in the 'palettes' directory).
            
            I shall go through both files with a fine toothed comb and explain
            as much as possible. (each method will be looked at carefully - if
            need be).            

==============================================================================
                        
Abstract view of the prototype:
-------------------------------

OK, start the XVT GUI, GUI then spawns the Geant benchmark program.  The user
interacts with the GUI via appropriate widgets selecting and executing
commands, these are sent to Geant which does the processing, Geant then 
relays any information back to the user via the GUI.

==============================================================================

More Detailed view:
-------------------

My main idea(s) for this interface:
            
Basically, the application starts (by executing the XVT GUI
executable file).
            
This handles all the initialisation required for the IPC mechanism,
it is in charge of creating the named pipes and opens them at the
XVT end with the appropriate permissions on them.
            
Then it handles the startup of the Geant 'xvt.benchmark' file.
            
It basically forks off a process and execs the xvt.benchmark file.
            
See the xvt.benchmark documentation for the ideas behind that and a
detailed description of what happens when it is executed.

==============================================================================

I will now cover some of the major parts and ideas of the interface before
analysing the individual methods:

==============================================================================

The named pipes:
----------------

This was the chosen IPC method after reviewing several IPC mechanisms (See 
the IPC section on my Geant4 web page).  It works like the twoWay 
program - incorporating a timer in the XVT end to read the pipe, a suitable 
way of testing the pipe before reading it (blocking issues) etc - see the 
early work in the IPC section or the program code for details.

The data Structures:
--------------------

Now, this was a new idea.  I had to get all the command information from Geant
and transfer it via the pipes and then use it to build a suitable dynamic
widget at the XVT end.  Therefore, the idea behind it all is simple, build a
data structure at the Geant end and fill it with command information, then when
requested by XVT, send all the information down the pipe.  At the XVT end, parse
this information into a data structure at this end and then build the widget
etc using this data structures contained information.

The new data structures I defined were G4CommandData and G4ParameterData. (See
the header file "g4xvtwin.h" for implementation specifics.

Some important issues were encountered.  It became clear that a certain 
'protocol' had to be defined so that I could recognise command information etc 
- this will be explained when documenting the individual method in detail.

The Palettes:
-------------

These were a new provided widget with XVT version 4.5 and were exactly what I
required. (I was working to emulate their behaviour before the release came
out).  The best feature of them is that they are dynamic and can be created at
run time etc - thus giving the flexibility needed for Geant.

As they are like an artists palette the first point is that you make a
selection and it stays selected - this suited my idea well as I had the idea
to have all the command information in one palette - i.e. select a command from
this palette, then have a button which you press to 'Execute' the selected
command.  I have since added a 'Help' button which you can also press - this
gives a small message about the selected command.

 - This is my idea for command entry to date:
 
 1.. Select command.
 2.. If help is required click help button for information.
 3.. If you are sure, press Execute which will dispatch the command to Geant.
 
Thus it gives a level of security so that selecting a command doesn't
necessarily mean it gets executed - until you EXPLICITLY say so.

Command Retrieval:
------------------

The idea behind this is quite clever.  OK, firstly, on the geant side of
things.  This reads all the commands from the Geant command tree using 'get'
methods from the Command tree class.  All of the commands are read into a type
I have defined - called G4CommandData - this basically holds all the 
information on a command and its parameters. (see data structures section
mentioned earlier)

Now back to the GUI end.  After this has initialised itself and listed all the
available commands the user can click on a button to retrieve all the available
commands from Geant.  When this happens a command is sent from the GUI down the
named pipes to Geant.  When Geant receives this command it parses the command
array at that end and writes all the information down the named pipe to XVT.

XVT also has the same data structure declared and basically reads all the
information from the pipe into its own structure.

After this is done the floating command palette is constructed using the
information contained in the data structure at the XVT end.

I have defined a protocol so that XVT recognises when the commands are coming
down the pipe.

As the timer ticks and calls the DoTimer method this checks to see if there is
anything present in the pipe - if so reads it and decides what to do with it.

When it receives the string "comTree" it knows that the commands are coming.

I basically bypass the timer for a while and then call some methods to deal 
with the incoming data.  The first string that comes through is an integer 
which tells me exactly how many command entries there are - then a loop can 
be entered to fill the array at the XVT end with that number of entries.

All the appropriate casting and type changes are handled correctly. (i.e.
converting from string to integer).

You have to be careful when reading the pipe that you check there is something
present before actually reading it - especially in this case when order of
information is critical - thats why I defined a method (getString) to handle 
correct reading of the named pipe. (it waits until thers is something present
before reading and returning the string).

==============================================================================

Now I shall look at the header file g4xvtwin.h and then I shall look at
the implementation file g4xvtwin.cxx in detail to explain what the methods 
do etc.

==============================================================================

g4xvtwin.h:
-----------

First thing to note in this file is the #defines which define XvtToGeant and 
GeantToXvt. These are names I will use in my code for the named pipes which 
I set up for communication between the XVT executable and the Geant4 executable. 

The next #define simply defines the file mode (i.e. appropriate UNIX 
permissions) to open the named pipes with.

After this come several #includes which include various functions I will 
use in the code.

After this I define 2 new types called 'G4parameterData' and 'G4commandData'. 
These are structures to hold information about commands.

After this are the usual method declarations and data member definitions 
etc. - These are commented.

******************************************************************************
*************************** IMPORTANT  *************************************** 
***                                                                        ***
*** The declaration of 'G4commandData commandArray[200];' in the header    ***
*** file basically declares how many commands can be stored in the array   ***
*** and then used to create the command palette.                           ***
*** IF during future use the GUI coredumps when creating the command       ***
*** palette it is most likely that there are greater that 200 commands     ***
*** - especially as commands are being added all the time.  All you need   ***
*** to do is increase this value and everything will work fine again.      ***
***                                                                        ***
******************************************************************************

==============================================================================

g4xvtwin.cxx:
-------------

I shall cover the methods which I have written or added code to - the ones 
automatically generated are self explanatory and have good commenting (see 
file for more details).

The Constructor:
----------------

Initialisation occurs here by calling several initialiser functions - see their
own documentation later on (next).

initIndex:
----------

This initialises the command base array indexing variables - very simple.

initStatusBar:
--------------

This sets up the status bars fields with some initial values.

initTimer:
----------

Starts the XVT timer (on every time interval (which the coder can alter) it
activates the DoTimer method).

initIPC:
--------

This sets the value of the pipe reading variable (i.e. how many bytes are 
read in a single read operation and then calls createNamedPipes and startSlave.

initFrame:
----------

Creates an attachment frame and allows attachment to all sides.

createNamedPipes:
-----------------

This method creates the named pipes required in the application.
  
It first creates the pipes and then opens them.  

The XvtToGeant pipe with write permission only, 
the GeantToXvt with read permission only.

The program uses the errorHandler method to report errors because the
program should terminate if IPC setup fails.

destroyNamedPipes:
------------------

This method closes and unlinks all the open named pipes.

startSlave:
-----------

This method forks and exec's the other process (xvt.benchmark - the Geant
executable).
 
This method also uses the errorHandler to report errors because if the 
fork or Exec fail the program should terminate.

The fork and exec are UNIX specific - for more details on their use see a UNIX
manual or man pages.

DoTimer:
--------

Every time the timer ticks this method checks to see if the file descriptor 
can be read without blocking (i.e. there is something in the pipe), if 
there is it grabs the data from the pipe, if you are not in the process of 
reading all the commands from Geant it then checks if the string is equal to 
"comTree" which indicates all the commands are coming from Geant and 
different action should be taken.  If it isn't comTree then it simply 
displays the string in the Main log area.

The different action to take is to call some other methods to do the reading of
commands and building of the command palette. (discussed later.)

Refer to code for the implementation specifics.

errorHandler:
-------------

A general purpose error handler.  It can be called by any method.  A string
is passed so that this method can post a dialog box to the user stating 
where the program fault occurs then quit the GUI.

The key point is that this method should be called when processing should no
longer continue i.e. the GUI should be killed.

writeData:
----------

Writes the passed string to the Command log window and the XvtToGeant 
pipe - i.e. will send the string to Geant (commands!) for processing.

initCommandPalette:
-------------------

Allocate storage space for a tool palette.

destroyCommandPalette:
----------------------

Free up storage space allocated for the command palette.

buildCommandPalettes:
---------------------

This method builds the palettes of commands available to the user. 
 
This basically parses the commandArray that was received from Geant and
builds up the command palette 'on-the-fly' from this information.  

It checks the flag of the command in the array to see whether it is
a 'directory', a command or a cascade. (this is denoted by the flag being a
0, 1 or 2 and this is set at the Geant end when it reads in the commands).

Once it has this information it can decide how to handle it.. If it is a
'directory' - i.e. a category of command it must have its own sub-palette
so this can be constructed.  If it is a command it must be a button in the
current sub-palette and if it is a cascade it must be a sub-palette of the
current palette - see code for implementation specifics.

destroyDragSource:
------------------

This method destroys the space allocated for the drag source.

setStatus:
----------

This method places the passed string into the status field of the status 
bar.

resetStatus:
------------

This method completely blanks the status field of the status bar.

retrieveCommand:
----------------

This method checks the command palette to see which button has been pressed 
then it will scan the command array to find the command that corresponds to
the users selection and returns the command string to be sent to Geant.

It also tests to see if the command that is selected for execution requires a
parameter and if so calls 'getParameters' to obtain this from the user.

initTextWindows:
----------------

This method simply puts a title in the Text windows to start things off.

fillLocalArray:
---------------

This method fills the command data structure at the XVT end with data that is
being received from the Geant end.  It uses the method 'getString' to grab the
strings from the pipe without error.

showHelp:
---------

This method is called when the user clicks on the 'help' button.

It looks which is the current selected tool and retrieves the guidance for
that tool - then puts this in the text widget and a pop up window.

getParameters:
--------------

This method retrieves parameters from the user via dialog boxes. - called by
retrieveCommand.

getTextCommand:
---------------

This method retrieves parameters from the user via dialog boxes.  This method
is very similar to 'getParameters' however this one is called if the user
wishes to type a command in instead of using the command palette widget.

executeCommand:
---------------

When the user clicks on the execute button this method is called.

It retrieves the command (with parameters if necessary) and writes it
down the named pipe to Geant to deal with.

getString:
----------

This method grabs the string that is located in the pipe and returns it 
to the caller - provided there is something present in the pipe.

activateButtons:
----------------

This method enables the 'exec' buttons on the toolbar.  It is called just
after the user has retrieved the commands and built the widget. 

It first enables the buttons and then refreshes the screen by redrawing 
them.

I have tried to limit the mistakes a user could make and this is a very good
way of doing just that - by disabling buttons/widgets etc that could cause
problems if pressed/clicked on at the wrong time.

==============================================================================

NOTE: Any code not mentioned here is self explanatory and well commented.

==============================================================================

Now I a couple of XVT specific discoveries to mention:

DISCOVERY 1:

There is some sort of problem with Dialog boxes and floating dockable palettes.

I had a palette which contained a button, when this button was pressed it 
called one of the pre-defined dialogs, this however caused a lock up concerning
the mouse - I couldn't get control back, the dialog seemed to keep it.
What happened exactly was when you clicked the button down it didn't pop up 
again thus it wasn't receiving the mouse up event.  When you did click the 
mouse again however it spawned the dialog - again and again, it sort of 
entered a loop repeatedly spawning the dialog. - Technical support said that
the dialog must have been inadvertently destroying the mouse handlers.

DISCOVERY 2:

You cannot have a toolbar with buttons in and a CAttachmentFrame which 
encomapasses the whole window because if you do, the CAttachmentFrame
obscures all of the mouse events from reaching the toolbar - thus none
of the buttons react.  

A solution to this problem is to have a CAttachmentFrame that doesn't
cover the entire window, just a part of it (this can be acheived with
a certain constructor - see code) - thus the events can then
get to the toolbar.

DISCOVERY 3:

For some reason you cannot pass a WINDOW variable between a view and a
document - when I tried it the program kept core dumping on this.

DISCOVERY 4:

The Select() system call has been dropped.  The reason is that it was failing
for some reason which I couldn't work out.  It kept returning TRUE all the
time when I used it to see if the named pipe was ready for reading and this
was not always the case.
  
The solution I came up with was to use a function from the ioctl library. 
The function I used is the following: 
  
int nbytes; 
ioctl(fd_XvtToGeant, FIONREAD, &nbytes);
  
This is much better than select because I only have a single file descriptor
to test and using this is far shorter and easier to code/understand than the
select stuff.

==============================================================================
S.Prior - August '97
