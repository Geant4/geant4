// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UItcsh.cc,v 1.2 2000-06-14 03:19:00 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef WIN32

#include <ctype.h>
#include "g4std/strstream"
#include "G4StateManager.hh"
#include "G4UIcommandStatus.hh"
#include "G4UItcsh.hh"

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
}

/////////////////////
G4UItcsh::~G4UItcsh()
/////////////////////
{
}
  
///////////////////////////
void G4UItcsh::MakePrompt()
///////////////////////////
{
  if(promptSetting.length()<=1) {
    promptString= promptSetting;
    return;
  }

  promptString="";
  G4int i;
  for(i=0; i<promptSetting.length()-1; i++){
    if(promptSetting[i]=='%'){
      switch (promptSetting[i+1]) {
      case 's':  // current application status
	{
	G4StateManager* statM= G4StateManager::GetStateManager();
	G4String stateStr= statM-> GetStateString(statM->GetCurrentState());
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
	char st[20];
	G4std::ostrstream os(st,20);
        os << currentHistoryNo << '\0';
	promptString.append(st);
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
  G4int i;
  for(i=cursorPosition-1; i<commandLine.length() ;i++) 
    G4cout << commandLine[i];
  for(i=cursorPosition-1; i<commandLine.length() ;i++)
    G4cout << AsciiBS;
  G4cout << G4std::flush;
    
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
    G4cout << AsciiBS << ' ' << AsciiBS << G4std::flush;
  } else { 
    G4cout << AsciiBS;
    G4int i;
    for(i=cursorPosition-2; i< commandLine.length()-1 ;i++){
      G4cout << commandLine[i+1];
    }
    G4cout << ' ';
    for(i=cursorPosition-2; i< commandLine.length() ;i++){
      G4cout << AsciiBS;
    }
    G4cout << G4std::flush;
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
  G4int i;
  for(i=cursorPosition-1; i< commandLine.length()-1 ;i++){
    G4cout << commandLine[i+1];
  }
  G4cout << ' ';
  for(i=cursorPosition-1; i< commandLine.length() ;i++){
    G4cout << AsciiBS;
  }
  G4cout << G4std::flush;

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
  for(i=1; i<=commandLine.length(); i++) G4cout << ' ';
  for(i=1; i<=commandLine.length(); i++) G4cout << AsciiBS;
  G4cout << G4std::flush;
  
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
  for(i=cursorPosition; i<=commandLine.length(); i++) G4cout << ' ';
  for(i=commandLine.length(); i>=cursorPosition; i--) G4cout << AsciiBS;
  G4cout << G4std::flush;

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

    G4cout << promptString << commandLine << G4std::flush;
    // reset cursur position
    for(G4int i=commandLine.length()+1; i>=cursorPosition+1; i--) 
      G4cout << AsciiBS << G4std::flush;
  }
}

//////////////////////////////
void G4UItcsh::ForwardCursor()
//////////////////////////////
{
  if(IsCursorLast()) return;

  G4cout << commandLine[cursorPosition-1] << G4std::flush;
  cursorPosition++;
}

///////////////////////////////
void G4UItcsh::BackwardCursor()
///////////////////////////////
{
  if(cursorPosition==1) return;

  cursorPosition--;
  G4cout << AsciiBS << G4std::flush;
}

//////////////////////////////
void G4UItcsh::MoveCursorTop()
//////////////////////////////
{
  for(G4int i=cursorPosition; i>1; i--){
    G4cout << AsciiBS;
  }
  G4cout << G4std::flush;
  cursorPosition=1;
}

//////////////////////////////
void G4UItcsh::MoveCursorEnd()
//////////////////////////////
{
  for(G4int i=cursorPosition-1; i<commandLine.length(); i++){
    G4cout << commandLine[i];
  }
  G4cout << G4std::flush;
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

    G4cout << commandLine << G4std::flush;
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

    G4cout << commandLine << G4std::flush;      
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
  if(jhead != G4String::npos) {
    input.remove(0, jhead);
    input= input.strip(G4String::leading);
  }

  // command tree of "user specified directory"
  G4String vpath= currentCommandDir;
  G4String vcmd;

  if( !input.empty() ) {
    G4int len= input.length();
    G4int indx=-1;
    for(G4int i=len-1; i>=0; i--) {
      if(input[i]=='/') {
        indx= i;
        break;
      }   
    }
    // get abs. path
    if(indx != -1) vpath= GetAbsCommandDirPath(input(0,indx+1));  
    if(!(indx==0  && len==1)) vcmd= input(indx+1,len-indx-1);  // care for "/"
  }

  // list matched dirs/commands
  ListCommand(vpath, vpath+vcmd);

  G4cout << promptString << commandLine << G4std::flush;
}

////////////////////////////////
void G4UItcsh::CompleteCommand()
////////////////////////////////
{
  // inputting string
  G4String input= G4String(commandLine).strip(G4String::leading);
  // target token is last token
  G4int jhead= input.last(' ');
  if(jhead != G4String::npos) {
    input.remove(0, jhead);
    input= input.strip(G4String::leading);
  }

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

  if(nMatch==0) return;  // no matched

  // display...
  input= commandLine;
  // target token is last token
  jhead= input.last(' ');
  if(jhead == G4String::npos) jhead=0;
  else jhead++;

  G4int jt= input.find_last_of('/');
  if(jt<jhead) jt=G4String::npos;

  if(jt==G4String::npos) jt= jhead;
  else jt++;

  G4String dspstr; 
  G4int i;
  for(i=jt; i<=input.length()-1; i++) dspstr+= G4String(AsciiBS); // cleanup
  for(i=jt; i<=input.length()-1; i++) dspstr+= G4String(' '); 
  for(i=jt; i<=input.length()-1; i++) dspstr+= G4String(AsciiBS); 

  dspstr+= stream;
  G4cout << dspstr << G4std::flush; 

  // command line string
  input.remove(jt);
  input+= stream;

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
	G4cout << G4endl;
	exit(0);
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
      if (cc == '[') {
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

///////////////////////////////////
G4String G4UItcsh::GetCommandLine()
///////////////////////////////////
{
  SetTermToInputMode();

  MakePrompt(); // update
  relativeHistoryIndex= 0;

  G4cout << promptString << G4std::flush;

  G4String newCommand= ReadLine();  // read line...
  // multi-line
  while( newCommand[newCommand.length()-1] == '_' ) {
    newCommand.remove(newCommand.length()-1);
    G4cout << G4endl;
    promptString= "? ";
    G4cout << promptString << G4std::flush;
    G4String newLine= ReadLine();
    newCommand.append(newLine);
  }

  // update history...
  G4bool isMeaningfull= FALSE; // check NULL command
  for (int i=0; i<newCommand.length(); i++) {
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
  for(int i=0; i<nmin; i++){
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

