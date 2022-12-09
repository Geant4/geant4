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
//

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4StateManager.hh"
#include "G4UIcommandStatus.hh"
#include "G4VUIshell.hh"
#include "G4UIArrayString.hh"

// terminal color string
static const G4String strESC(1,'\033');
static const G4String TermColorString[8] ={ 
  strESC+"[30m", strESC+"[31m", strESC+"[32m", strESC+"[33m", 
  strESC+"[34m", strESC+"[35m", strESC+"[36m", strESC+"[37m"
};

///////////////////////////////////////////////////////////////////
G4VUIshell::G4VUIshell(const G4String& prompt)
  : promptSetting(prompt), promptString(""), nColumn(80),  
    lsColorFlag(FALSE), directoryColor(BLACK), commandColor(BLACK),
    currentCommandDir("/")
///////////////////////////////////////////////////////////////////
{
}

/////////////////////////
G4VUIshell::~G4VUIshell()
/////////////////////////
{
}

////////////////////////////////////////////
void G4VUIshell::MakePrompt(const char* msg) 
////////////////////////////////////////////
{
  if(promptSetting.length()<=1) {
    promptString= promptSetting;
    return;
  }

  promptString="";
  G4int i;
  for(i=0; i<(G4int)promptSetting.length()-1; ++i){
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
      default:
	promptString += promptSetting[i];
        break;
      }           
    } else {
      promptString += promptSetting[i];
    }
  }

  // append last chaacter
  if(i == G4int(promptSetting.length())-1) 
    promptString += promptSetting[i];
}


////////////////////////////////
void G4VUIshell::ResetTerminal()
////////////////////////////////
{

}

// --------------------------------------------------------------------
//      G4command operations
// --------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////
G4UIcommandTree* G4VUIshell::GetCommandTree(const G4String& input) const
////////////////////////////////////////////////////////////////////////
{
  G4UImanager* UI= G4UImanager::GetUIpointer();

  G4UIcommandTree* cmdTree= UI-> GetTree();  // root tree

  G4String absPath = GetAbsCommandDirPath(G4StrUtil::strip_copy(input));

  // parsing absolute path ...
  if(absPath.length()==0) return NULL;
  if(absPath[G4int(absPath.length()-1)] != '/') return NULL; // error??
  if(absPath=="/") return cmdTree;

  for(std::size_t indx=1; indx<absPath.length()-1; ) {
    std::size_t jslash= absPath.find("/", indx);  // search index begin with "/" 
    if(jslash != G4String::npos) {
      if(cmdTree != NULL)
        cmdTree= cmdTree-> GetTree(G4String(absPath.substr(0,jslash+1)));
    }
    indx= jslash+1;
  }

  if(cmdTree == NULL) return NULL;
  else return cmdTree;
}

//////////////////////////////////////////////////////////////////////
G4String G4VUIshell::GetAbsCommandDirPath(const G4String& apath) const
//////////////////////////////////////////////////////////////////////
{
  if(apath.empty()) return apath;  // null string

  // if "apath" does not start with "/", 
  //   then it is treared as relative path
  G4String bpath= apath;
  if(apath[(std::size_t)0] != '/') bpath= currentCommandDir + apath;

  // parsing...
  G4String absPath= "/";
  for(std::size_t indx=1; indx<=bpath.length()-1; ) {
    std::size_t jslash= bpath.find("/", indx);  // search index begin with "/"
    if(indx == jslash) { // skip first '///'
      ++indx;
      continue;
    }
    if(jslash != G4String::npos) {
      if(bpath.substr(indx,jslash-indx) == ".."){  // directory up
        if(absPath == "/") {
          indx = jslash+1;
          continue;
        }
        if(absPath.length() >= 2) {
          absPath.erase(absPath.length()-1);  // remove last  "/"
          auto jpre= absPath.rfind('/');
          if(jpre != G4String::npos) absPath.erase(jpre+1);
        }
      } else if(bpath.substr(indx,jslash-indx) == "."){  // nothing to do
      } else { // add
        if( !(jslash==indx && bpath[(G4int)indx]=='/') ) // truncate "////"
          absPath+= bpath.substr(indx, jslash-indx+1);
          // better to be check directory existence. (it costs!)
      }
      indx= jslash+1;
    } else { // directory ONLY (ignore non-"/" terminated string)
      break;
    }
  }

  return  absPath;
}


////////////////////////////////////////////////////////////////////
G4String G4VUIshell::GetCommandPathTail(const G4String& apath) const
////////////////////////////////////////////////////////////////////
{   // xxx/xxx/zzz -> zzz, trancate /// -> /
  if(apath.empty()) return apath;

  G4int lstr= (G4int)apath.length();

  // for trancating "/"
  G4bool Qsla= FALSE;
  if(apath[lstr-1]=='/') Qsla= TRUE;

  // searching last '/' from tail
  G4int indx= -1;
  for(G4int i=lstr-1; i>=0; --i) {
    if(Qsla && apath[i]!='/') Qsla= FALSE; // break "/" flag!!
    if(apath[i]=='/' && !Qsla) {
      indx= i;
      break;
    } 
  }

  if(indx==-1) return apath;  // not found

  if(indx==0  && lstr==1) { // "/"
    G4String nullStr;
    return nullStr;
  } else {  
    //G4String newPath= apath(indx+1,lstr-indx-1); 
    G4String newPath= apath;
    newPath= newPath.substr(indx+1,lstr-indx-1);
    return newPath;
  }
}

// --------------------------------------------------------------------
//      shell commands
// --------------------------------------------------------------------
/////////////////////////////////////////////////////////////
void G4VUIshell::ListCommand(const G4String& dir, 
			     const G4String& candidate) const
/////////////////////////////////////////////////////////////
{
  // specified directpry
  G4String input = G4StrUtil::strip_copy(dir); 

  // command tree of "user specified directory"
  G4String vpath= currentCommandDir;
  G4String vcmd;

  G4int len= (G4int)input.length();
  if(! input.empty()) {
    G4int indx= -1;
    for(G4int i=len-1; i>=0; --i) { // search last '/'
      if(input[i]=='/') {
        indx= i;
        break;
      }   
    }
    // get abs. path
    if(indx != -1) vpath= GetAbsCommandDirPath(input.substr(0,indx+1));
    if(!(indx==0  && len==1)) vcmd= input.substr(indx+1,len-indx-1); // care for "/"
  }

  // check "vcmd" is directory?
  G4String inputpath= vpath+vcmd;
  if(! vcmd.empty()){
    G4String tmpstr= inputpath + "/";
    if(GetCommandTree(tmpstr) != NULL) {
      vpath= tmpstr;
      vcmd= "";
    }
  }
      
  // check "vpath" directory exists?
  G4UIcommandTree* atree= GetCommandTree(vpath);  
  if(atree == NULL) {
    G4cout << "<" << input << ">: No such directory" << G4endl;
    return;
  }

  // list matched directories/commands
  G4String stream;
  G4bool isMatch= FALSE;

  G4int Ndir= atree-> GetTreeEntry();
  G4int Ncmd= atree-> GetCommandEntry();
  if(Ndir==0 && Ncmd==0) return;  // no contents
  
  // directory ...
  for(G4int idir=1; idir<=Ndir; idir++) {
    if(idir==1 && lsColorFlag) stream+= TermColorString[directoryColor];
    G4String fpdir= atree-> GetTree(idir)-> GetPathName();
    // matching test
    if(candidate.empty()) { // list all
      if(vcmd=="" || fpdir==inputpath) {
	stream+= GetCommandPathTail(fpdir); stream+= "  ";
	isMatch= TRUE;
      }
    } else { // list only matched with candidate
      if( fpdir.find(candidate, 0) == 0) {
	stream+= GetCommandPathTail(fpdir); stream+= "  ";
      }
    }
  }
  
  // command ...
  for(G4int icmd=1; icmd<=Ncmd; icmd++){
    if(icmd==1 && lsColorFlag) stream+= TermColorString[commandColor];
    G4String fpcmd= atree-> GetPathName() +
             atree-> GetCommand(icmd) -> GetCommandName();
    // matching test
    if(candidate.empty()) { // list all
      if(vcmd=="" || fpcmd==inputpath) {
	stream+= GetCommandPathTail(fpcmd); stream+= "*  ";
	isMatch= TRUE;
      }
    } else {  // list only matched with candidate
      if( fpcmd.find(candidate, 0) == 0) {
	stream+= GetCommandPathTail(fpcmd); stream+= "*  ";
      }
    }
  }
  
  // waring : not matched
  if(!isMatch && candidate.empty()) 
    G4cout << "<" << input 
	   << ">: No such directory or command" << std::flush;

  // display
  G4UIArrayString arrayString(stream);
  arrayString.Show(nColumn);
}

