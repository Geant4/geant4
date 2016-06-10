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
// $Id: G4UItcsh.cc 66892 2013-01-17 10:57:59Z gunter $
//

#ifndef WIN32

#include "G4Types.hh"
#include "G4StateManager.hh"
#include "G4UIcommandStatus.hh"
#include "G4UItcsh.hh"
#include <ctype.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>

// ASCII character code
static const char AsciiCtrA = '\001';
static const char AsciiCtrB = '\002';
static const char AsciiCtrC = '\003';
static const char AsciiCtrD = '\004';
static const char AsciiCtrE = '\005';
static const char AsciiCtrF = '\006';
static const char AsciiCtrK = '\013';
static const char AsciiCtrL = '\014';
static const char AsciiCtrN = '\016';
static const char AsciiCtrP = '\020';
static const char AsciiCtrQ = '\021';
static const char AsciiCtrS = '\023';
static const char AsciiCtrZ = '\032';
static const char AsciiTAB   = '\011';
static const char AsciiBS    = '\010';
static const char AsciiDEL   = '\177';
static const char AsciiESC   = '\033';

static const int AsciiPrintableMin = 32;

// history file
static const G4String historyFileName= "/.g4_hist";

/////////////////////////////////////////////////////////
G4UItcsh::G4UItcsh(const G4String& prompt, G4int maxhist)
  : G4VUIshell(prompt),
    commandLine(""), cursorPosition(1),
    commandHistory(maxhist), maxHistory(maxhist),
    currentHistoryNo(1), relativeHistoryIndex(0)
/////////////////////////////////////////////////////////
{  
  // get current terminal mode
  tcgetattr(0, &tios);

  // read a shell history file
  const char* path = getenv("HOME");
  if( path == NULL ) return;

  G4String homedir= path;
  G4String fname= homedir + historyFileName;

  std::ifstream histfile;
  enum { BUFSIZE= 1024 }; char linebuf[BUFSIZE];

  histfile.open(fname, std::ios::in);
  while (histfile.good()) {
    if(histfile.eof()) break;

    histfile.getline(linebuf, BUFSIZE);
    G4String aline= linebuf;
    aline.strip(G4String::both);
    if(aline.size() !=  0) StoreHistory(linebuf);
  }
  histfile.close();
}

/////////////////////
G4UItcsh::~G4UItcsh()
/////////////////////
{
  // store a shell history
  const char* path = getenv("HOME");
  if( path == NULL ) return;

  G4String homedir= path;
  G4String fname= homedir + historyFileName;

  std::ofstream histfile;
  histfile.open(fname, std::ios::out);

  G4int n0hist= 1;
  if( currentHistoryNo > maxHistory ) n0hist= currentHistoryNo-maxHistory+1;

  for (G4int i=n0hist; i<= currentHistoryNo; i++) {
    histfile << RestoreHistory(i) << G4endl;
  }
  
  histfile.close();
}

//////////////////////////////////////////
void G4UItcsh::MakePrompt(const char* msg)
//////////////////////////////////////////
{
  if(promptSetting.length()<=1) {
    promptString= promptSetting;
    return;
  }

  promptString="";
  size_t i;
  for(i=0; i<promptSetting.length()-1; i++){
    if(promptSetting[i]=='%'){
      switch (promptSetting[i+1]) {
      case 's':  // current application status
	{
          G4String stateStr;
          if(msg)
          { stateStr = msg; }
          else
          {
	    G4StateManager* statM= G4StateManager::GetStateManager();
	    stateStr= statM-> GetStateString(statM->GetCurrentState());
          }
	  promptString.append(stateStr);
	  i++;
	}
        break;
      case '/':  // current working directory
	promptString.append(currentCommandDir);
	i++;
        break;
      case 'h':  // history#
	{
	std::ostringstream os;
        os << currentHistoryNo;
	promptString.append(os.str());
	i++;
	}
        break;
      default:
        break;
      } 
    } else {
      promptString.append(G4String(promptSetting[i]));
    }
  }

  // append last chaacter
  if(i == promptSetting.length()-1) 
    promptString.append(G4String(promptSetting[i]));
}


//////////////////////////////
void G4UItcsh::ResetTerminal()
//////////////////////////////
{
  RestoreTerm();
}


// --------------------------------------------------------------------
//      commad line operations
// --------------------------------------------------------------------
//////////////////////////////////////
void G4UItcsh::InitializeCommandLine()
//////////////////////////////////////
{
  commandLine= "";
  cursorPosition= 1;
}

///////////////////////////////////////
void G4UItcsh::InsertCharacter(char cc)
///////////////////////////////////////
{
  if( ! (cc >= AsciiPrintableMin  && isprint(cc)) ) return;

  // display...
  G4cout << cc;
  size_t i;
  for(i=cursorPosition-1; i<commandLine.length() ;i++) 
    G4cout << commandLine[i];
  for(i=cursorPosition-1; i<commandLine.length() ;i++)
    G4cout << AsciiBS;
  G4cout << std::flush;
    
  // command line string...
  if(IsCursorLast()) {  // add
    commandLine+= cc;
  } else { // insert
    commandLine.insert(cursorPosition-1, G4String(cc));
  }
  cursorPosition++;
}
  
///////////////////////////////////
void G4UItcsh::BackspaceCharacter()
///////////////////////////////////
{
  if(cursorPosition==1) return;

  // display...
  if(IsCursorLast()) {  
    G4cout << AsciiBS << ' ' << AsciiBS << std::flush;
  } else { 
    G4cout << AsciiBS;
    size_t i;
    for(i=cursorPosition-2; i< commandLine.length()-1 ;i++){
      G4cout << commandLine[i+1];
    }
    G4cout << ' ';
    for(i=cursorPosition-2; i< commandLine.length() ;i++){
      G4cout << AsciiBS;
    }
    G4cout << std::flush;
  }

  // command line string...
  commandLine.erase(cursorPosition-2, 1);

  cursorPosition--;
}

////////////////////////////////
void G4UItcsh::DeleteCharacter()
////////////////////////////////
{
  if(IsCursorLast()) return;

  // display...
  size_t i;
  for(i=cursorPosition-1; i< commandLine.length()-1 ;i++){
    G4cout << commandLine[i+1];
  }
  G4cout << ' ';
  for(i=cursorPosition-1; i< commandLine.length() ;i++){
    G4cout << AsciiBS;
  }
  G4cout << std::flush;

  // command lin string...
  commandLine.erase(cursorPosition-1, 1);
}

//////////////////////////
void G4UItcsh::ClearLine()
//////////////////////////
{
  // display...
  G4int i;
  for(i= cursorPosition; i>=2; i--) G4cout << AsciiBS;
  for(i=1; i<=G4int(commandLine.length()); i++) G4cout << ' ';
  for(i=1; i<=G4int(commandLine.length()); i++) G4cout << AsciiBS;
  G4cout << std::flush;
  
  // command line string...
  commandLine.erase();
  cursorPosition= 1;
}

/////////////////////////////////
void G4UItcsh::ClearAfterCursor()
/////////////////////////////////
{
  if(IsCursorLast()) return;

  // display...
  G4int i;
  for(i=cursorPosition; i<=G4int(commandLine.length()); i++) G4cout << ' ';
  for(i=commandLine.length(); i>=cursorPosition; i--) G4cout << AsciiBS;
  G4cout << std::flush;

  // command line string...
  commandLine.erase(cursorPosition-1, 
		    commandLine.length()-cursorPosition+1);
}

////////////////////////////
void G4UItcsh::ClearScreen()
////////////////////////////
{
  if(! clearString.empty() ) {
    G4cout << clearString;

    G4cout << promptString << commandLine << std::flush;
    // reset cursur position
    for(G4int i=commandLine.length()+1; i>=cursorPosition+1; i--) 
      G4cout << AsciiBS << std::flush;
  }
}

//////////////////////////////
void G4UItcsh::ForwardCursor()
//////////////////////////////
{
  if(IsCursorLast()) return;

  G4cout << commandLine[(size_t)(cursorPosition-1)] << std::flush;
  cursorPosition++;
}

///////////////////////////////
void G4UItcsh::BackwardCursor()
///////////////////////////////
{
  if(cursorPosition==1) return;

  cursorPosition--;
  G4cout << AsciiBS << std::flush;
}

//////////////////////////////
void G4UItcsh::MoveCursorTop()
//////////////////////////////
{
  for(G4int i=cursorPosition; i>1; i--){
    G4cout << AsciiBS;
  }
  G4cout << std::flush;
  cursorPosition=1;
}

//////////////////////////////
void G4UItcsh::MoveCursorEnd()
//////////////////////////////
{
  for(size_t i=cursorPosition-1; i<commandLine.length(); i++){
    G4cout << commandLine[i];
  }
  G4cout << std::flush;
  cursorPosition=commandLine.length()+1;
}

////////////////////////////////
void G4UItcsh::PreviousCommand()
////////////////////////////////
{
  G4int nhmax= currentHistoryNo-1 >= maxHistory ? 
                 maxHistory : currentHistoryNo-1;

  // retain current input
  if(relativeHistoryIndex==0) commandLineBuf= commandLine;

  if(relativeHistoryIndex>=-nhmax+1 && relativeHistoryIndex<=0) {
    ClearLine();
    relativeHistoryIndex--;
    commandLine= RestoreHistory(currentHistoryNo+relativeHistoryIndex);

    G4cout << commandLine << std::flush;
    cursorPosition= commandLine.length()+1;
  }
}

////////////////////////////
void G4UItcsh::NextCommand()
////////////////////////////
{  
  G4int nhmax= currentHistoryNo-1 >= maxHistory ? 
                 maxHistory : currentHistoryNo-1;

  if(relativeHistoryIndex>=-nhmax && relativeHistoryIndex<=-1) {
    ClearLine();
    relativeHistoryIndex++;

    if(relativeHistoryIndex==0) commandLine= commandLineBuf;
    else commandLine= RestoreHistory(currentHistoryNo+relativeHistoryIndex);

    G4cout << commandLine << std::flush;      
    cursorPosition= commandLine.length()+1;
  }
}


///////////////////////////////////
void G4UItcsh::ListMatchedCommand()
///////////////////////////////////
{
  G4cout << G4endl;
  
  // input string
  G4String input= G4String(commandLine).strip(G4String::leading);
  // target token is last token
  G4int jhead= input.last(' ');
  if(jhead != G4int(G4String::npos)) {
    input.remove(0, jhead);
    input= input.strip(G4String::leading);
  }
  //G4cout << "@@@@ input=" << input << G4endl;

  // command tree of "user specified directory"
  G4String vpath = currentCommandDir;
  G4String vcmd = "";

  if( !input.empty() ) {
    G4int len= input.length();
    G4int indx=-1;
    for(G4int i=len-1; i>=0; i--) {
      if(input[(size_t)i]=='/') {
        indx= i;
        break;
      }
    }
    // get abs. path
    if(indx != -1) vpath= GetAbsCommandDirPath(input(0,indx+1));  
    if(!(indx==0  && len==1)) vcmd= input(indx+1,len-indx-1);  // care for "/"
  }

  // list matched dirs/commands
  //G4cout << "@@@ vpath=" << vpath <<":vcmd=" << vcmd << G4endl;
  ListCommand(vpath, vpath+vcmd);

  G4cout << promptString << commandLine << std::flush;
}

////////////////////////////////
void G4UItcsh::CompleteCommand()
////////////////////////////////
{
  // inputting string
  G4String input= G4String(commandLine).strip(G4String::leading);
  // target token is last token
  G4int jhead= input.last(' ');
  if(jhead != G4int(G4String::npos)) {
    input.remove(0, jhead);
    input= input.strip(G4String::leading);
  }

  // tail string
  size_t thead = input.find_last_of('/');
  G4String strtail = input;
  if (thead != G4String::npos) strtail = input(thead+1, input.size()-thead-1);

  // command tree of "user specified directory"  
  G4String vpath= currentCommandDir;
  G4String vcmd;

  G4int len= input.length();
  if(!input.empty()) {
    G4int indx= -1;
    for(G4int i=len-1; i>=0; i--) {
      if(input(i)=='/') {
        indx= i;
        break;
      }   
    }
    // get abs. path
    if(indx != -1) vpath= GetAbsCommandDirPath(input(0,indx+1));
    if(!(indx==0  && len==1)) vcmd= input(indx+1,len-indx-1);  // care for "/"
  }

  G4UIcommandTree* atree= GetCommandTree(vpath);  // get command tree
  if(atree == NULL) return;

  // list matched directories/commands
  G4String stream, strtmp;
  G4String inputpath= vpath+vcmd;
  G4int nMatch= 0;

  int Ndir= atree-> GetTreeEntry();
  int Ncmd= atree-> GetCommandEntry();
  
  // directory ...
  for(G4int idir=1; idir<=Ndir; idir++) {
    G4String fpdir= atree-> GetTree(idir)-> GetPathName();
    // matching test
    if( fpdir.index(inputpath, 0) == 0) {
      if(nMatch==0) {
        stream= GetCommandPathTail(fpdir);
      } else {
        strtmp= GetCommandPathTail(fpdir);
        stream= GetFirstMatchedString(stream, strtmp);
      }
      nMatch++;
    }
  }
  
  // command ...
  for(G4int icmd=1; icmd<=Ncmd; icmd++){
    G4String fpcmd= atree-> GetPathName() +
                    atree-> GetCommand(icmd) -> GetCommandName();
    // matching test
    if( fpcmd.index(inputpath, 0) ==0) {
      if(nMatch==0) {
        stream= GetCommandPathTail(fpcmd) + " ";
      } else {
        strtmp= GetCommandPathTail(fpcmd) + " ";
        stream= GetFirstMatchedString(stream, strtmp);
      }
      nMatch++;
    }
  }

  // display...
  input= commandLine;
  // target token is last token
  jhead= input.last(' ');
  if(jhead == G4int(G4String::npos)) jhead=0;
  else jhead++;

  G4int jt = jhead;

  G4String dspstr; 
  G4int i;
  for(i=jt; i<=G4int(input.length())-1; i++) dspstr+= G4String(AsciiBS); 
  for(i=jt; i<=G4int(input.length())-1; i++) dspstr+= G4String(' '); 
  for(i=jt; i<=G4int(input.length())-1; i++) dspstr+= G4String(AsciiBS); 

  dspstr+= (vpath + stream);
  if (nMatch == 0) dspstr+= strtail;
  G4cout << dspstr << std::flush;

  // command line string
  input.remove(jt);
  input+= (vpath + stream);
  if (nMatch==0) input+= strtail;

  commandLine= input;
  cursorPosition= commandLine.length()+1;
}

// --------------------------------------------------------------------
//      commad line
// --------------------------------------------------------------------
/////////////////////////////
G4String G4UItcsh::ReadLine()
/////////////////////////////
{
  InitializeCommandLine();

  char cc;
  do{  // input loop
    G4cin.get(cc);

    // treatment for special character
    switch(cc){
    case AsciiCtrA:       // ... move cursor to the top
      MoveCursorTop();
      break;
    case AsciiCtrB:       // ... backward cursor
      BackwardCursor();
      break;
    case AsciiCtrD:       // ... delete/exit/show matched list
      if(commandLine.length()!=0 && IsCursorLast()) ListMatchedCommand();
      else if (commandLine.empty()) {
        return G4String("exit");	
      } else DeleteCharacter();
      break;
    case AsciiCtrE:       // ... move cursor to the end
      MoveCursorEnd();
      break;
    case AsciiCtrF:       // ... forward cursor
      ForwardCursor();
      break;
    case AsciiCtrK:       // ... clear after the cursor
      ClearAfterCursor();
      break;
    case AsciiCtrL:       // ... clear screen
      // ClearScreen();
      break;
    case AsciiCtrN:	// ... next command
      NextCommand();
      break;
    case AsciiCtrP:	// ... previous command
      PreviousCommand();
      break;
    case AsciiTAB:         // ... command completion
      if( (!commandLine.empty()) && IsCursorLast()) CompleteCommand();
      break;
    case AsciiDEL:         // ... backspace
      BackspaceCharacter();
      break;
    case AsciiBS:          // ... backspace
      BackspaceCharacter();
      break;
    case AsciiCtrC:       // ... kill prompt
      break;
    case AsciiCtrQ:       // ... restarts suspeded output
      break;
    case AsciiCtrS:       // ... suspend output
      break;
    case AsciiCtrZ:       // ... suspend
      break;
    default:
      break;
    }

    // treatment for ESC. character
    if( cc == AsciiESC) { // ESC
      G4cin.get(cc);
      if (cc == '[' || cc == 'O') { // care for another termcap, such as konsole
        G4cin.get(cc);
        switch(cc) {
        case 'A': // [UP]
          cc = 'P' - '@';
          PreviousCommand();  // ... show previous commad
          break;
        case 'B': // [DOWN]
          cc = 'N' - '@';
          NextCommand();  // ... show next commad
          break;
        case 'C': // [RIGHT]
          cc = 'F' - '@';
          ForwardCursor();   // ... forward cursor
          break;
        case 'D': // [LEFT]
          cc = 'B' - '@';
          BackwardCursor();      // ... backward cursor
          break;
        default:  // who knows !?
          cc = 0;
          break;
        }
      }
    }

    // insert character to command line and display
    InsertCharacter(cc);
  
  } while( cc != '\n');

  return commandLine;
}

////////////////////////////////////////////////////////
G4String G4UItcsh::GetCommandLineString(const char* msg)
////////////////////////////////////////////////////////
{
  SetTermToInputMode();

  MakePrompt(msg); // update
  relativeHistoryIndex= 0;

  G4cout << promptString << std::flush;

  G4String newCommand= ReadLine();  // read line...
  // multi-line
  while( (newCommand.length() > 0) &&
	 ( newCommand[newCommand.length()-1] == '_') ) {
    newCommand.remove(newCommand.length()-1);
    G4cout << G4endl;
    promptString= "? ";
    G4cout << promptString << std::flush;
    G4String newLine= ReadLine();
    newCommand.append(newLine);
  }

  // update history...
  G4bool isMeaningfull= FALSE; // check NULL command
  for (size_t i=0; i<newCommand.length(); i++) {
    if(newCommand[i] != ' ') {
      isMeaningfull= TRUE;
      break;
    }
  }
  if( !newCommand.empty() && isMeaningfull) StoreHistory(newCommand);

  // reset terminal
  RestoreTerm();

  G4cout << G4endl;
  return newCommand;
}

////////////////////////////////////////////////////////////////////
G4String G4UItcsh::GetFirstMatchedString(const G4String& str1, 
					 const G4String& str2) const
////////////////////////////////////////////////////////////////////
{
  int nlen1= str1.length();
  int nlen2= str2.length();

  int nmin = nlen1<nlen2 ? nlen1 : nlen2;

  G4String strMatched;
  for(size_t i=0; G4int(i)<nmin; i++){
    if(str1[i]==str2[i]) {
      strMatched+= str1[i];
    } else {
      break;
    }
  }

  return strMatched;
}

// --------------------------------------------------------------------
//      history
// --------------------------------------------------------------------
//////////////////////////////////////////////
void G4UItcsh::StoreHistory(G4String aCommand)
//////////////////////////////////////////////
{
  G4int i= currentHistoryNo%maxHistory; 
  if(i==0) i=maxHistory;

  commandHistory[i-1]= aCommand;  // 0-offset
  currentHistoryNo++;
}

///////////////////////////////////////////////
G4String G4UItcsh::RestoreHistory(G4int histNo)
///////////////////////////////////////////////
{
  if(histNo>= currentHistoryNo) return "";

  G4int index= histNo%maxHistory;
  if(index==0) index= maxHistory;

  return commandHistory[index-1]; // 0-offset
}

// --------------------------------------------------------------------
//      terminal mode
// --------------------------------------------------------------------
///////////////////////////////////
void G4UItcsh::SetTermToInputMode()
///////////////////////////////////
{  
  termios tiosbuf= tios;

  tiosbuf.c_iflag &= ~(BRKINT | ISTRIP);
  tiosbuf.c_iflag |= (IGNBRK | IGNPAR);
  tiosbuf.c_lflag &= ~(ICANON | IEXTEN | ECHO);
  tiosbuf.c_cc[VMIN] = 1;
  tiosbuf.c_cc[VTIME] = 0;
  
  tcsetattr(0, TCSAFLUSH, &tiosbuf);
}


////////////////////////////
void G4UItcsh::RestoreTerm()
////////////////////////////
{
  tcsetattr(0, TCSAFLUSH, &tios);
}  

#endif

