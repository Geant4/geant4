// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: g4xvtwin.h,v 1.1 1999-01-07 16:05:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/********************************************************************
 *
    Copyright (c) 1996

    File Name:    G4XvtWin.h
    Author:       Simon Prior
    Date:

    Class:        G4XvtWin
    Inheritance:  CWindow->G4XvtWin
    Helper(s):

    Purpose:      This is a window class that has properties that were
                  specified using XVT-Architect. It implements methods
                  to Initialize and manage the state of views nested
                  inside the window.

    Usage:        Override or add methods as necessary.

    Modifications:
 *
 ********************************************************************/

#ifndef G4XvtWin_H
#define G4XvtWin_H

////////////////////////////////////////////////////////////////////////
// Define named pipes                                                 //
//                                                                    //
// Note - if you wish to change the location of where the pipes are   //
//        created edit the following lines putting in the appropriate //
//        paths.  The default is for the pipes to be created in the   //
//        current directory.                                          //
////////////////////////////////////////////////////////////////////////

#define XvtToGeant "XvtToGeant.tmp"  
#define GeantToXvt "GeantToXvt.tmp"


////////////////////////////////////////////////////////////////////////
// Define file mode/permissions                                       //
////////////////////////////////////////////////////////////////////////

#define FILE_MODE (0664 | S_IFIFO)   


////////////////////////////////////////////////////////////////////////
// Include different files depending on what platform you are         //
// compiling on.                                                      //
////////////////////////////////////////////////////////////////////////

#ifdef _AIX              
#include <sys/mode.h>    // --- mkfifo utils --- //
#endif

#include <strings.h>     // --- bzero function --- //
#include <sys/time.h>    // --- struct timeval --- //
#include <sys/types.h>
#include <unistd.h>      // --- POSIX standards for stdin, out, err. --- //
#include <values.h>
#include <sys/stat.h>
#include <fcntl.h>       // --- Open, Close, Unlink etc --- //
#include <sys/ioctl.h>
#include <stdio.h>       // --- sscanf, sprintf --- //

#include "XVTPwr.h"
#include "AppDef.h"

#include GUI_i
#include CWindow_i

#ifdef DSP_RELEASE
#include G4XvtWin_f
#endif


////////////////////////////////////////////////////////////////////////
// New data structures, one to hold the parameter information of a    //
// command and one to hold information on the command itself.         //
////////////////////////////////////////////////////////////////////////

typedef struct G4parameterData
{
  CStringRW name;
  CStringRW guidance;
  CStringRW defaultValue;
  CStringRW range;
  CStringRW candidate;
  char      type;
  CStringRW omittable; 
};

typedef struct G4commandData
{
  int flag;
  CStringRW name;
  CStringRW guidance;
  int numOfParameters;
  G4parameterData parameters[10];  

//  CStringRW bitmapName;

  int commandNumber;
};


////////////////////////////////////////////////////////////////////////
// Variables for the status bar                                       //
////////////////////////////////////////////////////////////////////////

const long NULL_FIELD = 1L;
const long STATUS_FIELD = 2L;
const long COORDS_FIELD = 3L;


////////////////////////////////////////////////////////////////////////
// Make code aware of the attachment and mouse handling classes use   //
////////////////////////////////////////////////////////////////////////

class CAttachmentFrame;
class CAttachment;
class CToolBarAttachment;
class CToolPalette;
class CViewSinkData;
class CDragSource;

class G4XvtWin : public CWindow
{
public:
	G4XvtWin(
		CDocument *      theDocument,
		const CRect&     theRegion,
		const CStringRW& theTitle            = NULLString,
		long             theWindowAttributes = WSF_NONE,
		WIN_TYPE         theWindowType       = W_DOC,
		int              theMenuBarId        = MENU_BAR_RID,
		WINDOW           theParentWindow     = TASK_WIN );

	virtual ~G4XvtWin( void );

	BOOLEAN  I4XvtWin( void );

	virtual void DoCommand( long theCommand, void * theData = NULL );

	virtual void DoMenuCommand( MENU_TAG theMenuItem,
		BOOLEAN isShiftKey, BOOLEAN isControlKey );

	virtual void UpdateMenus( CMenuBar * theMenuBar);

	virtual void DoUpdateModel( long theControllerId,
		long theCommand, const CModel * theModel );

        ////////////////////////////////////////////////////////////////
        // Overriden timer method                                     //
        ////////////////////////////////////////////////////////////////
    
        void DoTimer(long theTimerId);

protected:
	// By default, copy and assignment are disallowed

	G4XvtWin( const G4XvtWin& theWindow ) : CWindow( theWindow ) { }
	G4XvtWin& operator=( const G4XvtWin& theWindow ){ return( *this ); }
                        
private:
	G4WinData itsData;

        ////////////////////////////////////////////////////////////////
        // IPC Data members                                           //
        ////////////////////////////////////////////////////////////////
        
        int fd_XvtToGeant, fd_GeantToXvt;
        long theTimer;
        long interval;
        WINDOW theWin;
        char textDump[100];
        int number;   
        int G4Running;
        
        ///////////////////////////////////////////////////////////////       
        // IPC methods                                               //
        ///////////////////////////////////////////////////////////////
        
        void initTimer(void);
        void initIPC(void);
        void startSlave(void);
        void createNamedPipes(void);  
        void destroyNamedPipes(void);
        void errorHandler(CStringRW theError);
        void writeData(CStringRW theString);
        void ExecuteCommand(void);
        CStringRW getString(void);
        
        ///////////////////////////////////////////////////////////////       
        // Command Base data members                                 //
        ///////////////////////////////////////////////////////////////
        
        G4commandData commandArray[200];
        int currentPosition, lastPosition;
        int reading;
        
        ///////////////////////////////////////////////////////////////       
        // Command Base methods                                      //
        ///////////////////////////////////////////////////////////////
                
        void initIndex(void);
        void fillLocalArray(void);
        void listCommandArray(void); // debugging method       
        
        
        ///////////////////////////////////////////////////////////////       
        // Palette and attachment data members                       //
        ///////////////////////////////////////////////////////////////
                
        CToolPalette* commandPalette;
        CDragSource* commandSource;
        CAttachmentFrame* itsAttachmentFrame;
       
        ///////////////////////////////////////////////////////////////       
        // Palette and attachment methods                            //
        ///////////////////////////////////////////////////////////////
        
        void initFrame(void);
        void initCommandPalette(void);
        void destroyCommandPalette(void);
        void destroyDragSource(void);
        void buildCommandPalettes(void);                
        void regenerateCommandPalette(void);
        CStringRW retrieveCommand(void);

        
        ///////////////////////////////////////////////////////////////       
        // Status bar methods                                        //
        ///////////////////////////////////////////////////////////////
        
        void initStatusBar(void);       
        void setStatus(CStringRW theText);
        void resetStatus(void);
                
        
        ///////////////////////////////////////////////////////////////       
        // Text window methods                                       //
        ///////////////////////////////////////////////////////////////
        
        void initTextWindows(void);


        ///////////////////////////////////////////////////////////////       
        // Help Method                                               //
        ///////////////////////////////////////////////////////////////

        void showHelp(void);
 
        
        ///////////////////////////////////////////////////////////////       
        // Parameter retrieval methods (pops up dialog)              //
        ///////////////////////////////////////////////////////////////

        CStringRW getParameters(void);
        CStringRW getTextCommand(void);

        ///////////////////////////////////////////////////////////////       
        // Toolbar and enclosed view methods                         //
        ///////////////////////////////////////////////////////////////
        
        void activateButtons(void);

//    PWRClassInfo
};

#endif // G4XvtWin_H

