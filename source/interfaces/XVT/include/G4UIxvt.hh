// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIxvt.hh,v 1.2 1999-04-13 01:24:58 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
//                     XVT driver class header                        //
//                     ~~~~~~~~~~~~~~~~~~~~~~~                        //
// Written by: Simon Prior                                            //
//       Date: 22/04/97                                               //
//                                                                    //
// Updated for state machine: 12/08/97                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////
#ifndef G4UIxvt_h
#define G4UIxvt_h

#if defined(G4UI_BUILD_XVT_SESSION) || defined(G4UI_USE_XVT)

////////////////////////////////////////////////////////////////////////
// define the named pipe names                                        //
////////////////////////////////////////////////////////////////////////

#define XvtToGeant "XvtToGeant.tmp"  
#define GeantToXvt "GeantToXvt.tmp"


////////////////////////////////////////////////////////////////////////
// mode for opening a named pipe                                      //
////////////////////////////////////////////////////////////////////////

#define FILE_MODE (0664 | S_IFIFO)    


////////////////////////////////////////////////////////////////////////
// If compiling using AFS different libraries required                //
////////////////////////////////////////////////////////////////////////

#ifdef _AIX              
#include <sys/mode.h>    // mkfifo utils 
#endif

#include <unistd.h>      // POSIX standards for stdin, out, err.
//#include <values.h>
#include <sys/stat.h>    // mkfifo utils for HP 
#include <fcntl.h>       // Open, Close, Unlink etc 
#include <sys/ioctl.h>
#include <stdio.h>       // sscanf, sprinf 

#ifndef FIONREAD
#include <sys/filio.h>   // FIONREAD defines
#endif

////////////////////////////////////////////////////////////////////////
// Geant4 specific includes                                           //
////////////////////////////////////////////////////////////////////////

#include "G4UIsession.hh"
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"


////////////////////////////////////////////////////////////////////////
// State machine Geant4 includes                                      //
////////////////////////////////////////////////////////////////////////

#include "G4VStateDependent.hh"
#include <fstream.h>


////////////////////////////////////////////////////////////////////////
// New data structures, one to hold the parameter information of a    //
// command and one to hold information on the command itself.         //
////////////////////////////////////////////////////////////////////////

typedef struct G4parameterData
{
  G4String name;
  G4String guidance;
  G4String defaultValue;
  G4String range;
  G4String candidate;
  char     type;
  G4String omittable; 
} G4parameterData;

typedef struct G4commandData
{
  int flag;
  G4String name;
  G4String guidance;
  int numOfParameters;
  G4parameterData parameters[10];
  
//  G4String bitmapName;

  int commandNumber;
} G4commandData;


class G4UIxvt : public G4UIsession//, public G4VStateDependent
{
  public:
    G4UIxvt();
    ~G4UIxvt();
    
    ////////////////////////////////////////////////////////////////////
    // Inherited virtual functions                                    //
    ////////////////////////////////////////////////////////////////////

//    void SessionStart(void);    
G4UIsession*  SessionStart(void);    
    G4String GetCommand(void);
    void SessionTerminate(void);


    ////////////////////////////////////////////////////////////////////
    // New state machine methods                                      //
    ////////////////////////////////////////////////////////////////////
    
    G4bool Notify(G4ApplicationState requestedState);

    inline void set_breakPointAt(int id,G4bool flg)
    { 
      brktbl[id] = flg; 
      noBreakFlag = !(brktbl[0]||brktbl[1]||brktbl[2]||brktbl[3]);
    };

    void set_verbose(void);    
    
  private:
  
    ////////////////////////////////////////////////////////////////////
    // Terminal and interface data members (might be taken out as     //   
    // they are not needed for the XVT GUI)                           //
    ////////////////////////////////////////////////////////////////////    
        
    G4UImanager * UI;          


    ////////////////////////////////////////////////////////////////////
    // New state machine data members                                 //
    ////////////////////////////////////////////////////////////////////
    
    G4bool iExit;
    G4bool iCont;
    G4bool brktbl[4];
    G4bool noBreakFlag;
    G4bool verboseFlag;


    ////////////////////////////////////////////////////////////////////
    // New state machine methods                                      //
    ////////////////////////////////////////////////////////////////////
    
    void additionalSession(void);    
    G4bool breakRequested(G4ApplicationState,G4ApplicationState);


    ////////////////////////////////////////////////////////////////////
    // IPC Data Members                                               //
    ////////////////////////////////////////////////////////////////////
    
    int fd_XvtToGeant, fd_GeantToXvt; 
    char textDump[100];  
    int number;   


    ////////////////////////////////////////////////////////////////////
    // IPC methods                                                    //
    ////////////////////////////////////////////////////////////////////
    
    void openConnections(void);
    int checkXvtToGeantPipe(void);
    G4String readXvtToGeant(void);
    void writeGeantToXvt(const G4String& theString); 


    ////////////////////////////////////////////////////////////////////
    // Terminal and TCL methods updated for use with XVT              //
    ////////////////////////////////////////////////////////////////////
    
    void codeGen(G4UIcommandTree *,int recursive_level);                           
    void listCurrentDirectory(void);

    
    ////////////////////////////////////////////////////////////////////
    // Command handling data members                                  //
    ////////////////////////////////////////////////////////////////////
    
    G4commandData commandArray[200];
    int currentPosition, lastPosition;


    ////////////////////////////////////////////////////////////////////
    // Command handling methods                                       //
    ////////////////////////////////////////////////////////////////////
    
    void fillArrayEntries(G4UIcommandTree * tree, int level);
    void sendArrayEntriesToXvt(void);
    void briefListCommands(void);

        
    ////////////////////////////////////////////////////////////////////
    // Debugging methods                                              //
    ////////////////////////////////////////////////////////////////////
    void listCommandArray(void);
    void errorHandler(const G4String& theError);


};

#endif

#endif
