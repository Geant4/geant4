//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4UImanager.cc,v 1.15 2001-10-11 01:38:00 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------

#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIsession.hh"
#include "G4UIbatch.hh"
#include "G4UIcontrolMessenger.hh"
#include "G4UnitsMessenger.hh"
#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include "G4StateManager.hh"
#include "G4UIaliasList.hh"

#include "g4std/strstream"
#include "g4rw/ctoken.h"


G4UImanager * G4UImanager::fUImanager = 0;
G4bool G4UImanager::fUImanagerHasBeenKilled = false;

G4UImanager * G4UImanager::GetUIpointer()
{
  if(!fUImanager) 
  {
    if(!fUImanagerHasBeenKilled)
    {
      fUImanager = new G4UImanager;
      fUImanager->CreateMessenger();
    }
  }
  return fUImanager;
}

G4UImanager::G4UImanager()
:G4VStateDependent(true)
{
  savedCommand = 0;
  treeTop = new G4UIcommandTree("/");
  aliasList = new G4UIaliasList;
  G4String nullString;
  savedParameters = nullString;
  verboseLevel = 0;
  saveHistory = false;
  session = NULL;
  SetCoutDestination(session);
  pauseAtBeginOfEvent = false;
  pauseAtEndOfEvent = false;
}

void G4UImanager::CreateMessenger()
{
  UImessenger = new G4UIcontrolMessenger;
  UnitsMessenger = new G4UnitsMessenger;
}

G4UImanager::~G4UImanager()
{
  SetCoutDestination(NULL);
  histVec.clear();
  if(saveHistory) historyFile.close();
  delete UImessenger;
  delete UnitsMessenger;
  delete treeTop;
  delete aliasList;
  fUImanagerHasBeenKilled = true;
  fUImanager = NULL;
}

G4UImanager::G4UImanager(const G4UImanager &right) { }
const G4UImanager & G4UImanager::operator=(const G4UImanager &right)
{ return right; }
int G4UImanager::operator==(const G4UImanager &right) const
{ return false; }
int G4UImanager::operator!=(const G4UImanager &right) const
{ return true; }

G4String G4UImanager::GetCurrentValues(const char * aCommand)
{
  G4String theCommand = aCommand;
  savedCommand = treeTop->FindPath( theCommand );
  if( savedCommand == NULL ) 
  {
    G4cerr << "command not found" << G4endl;
    return G4String();
  }
  return savedCommand->GetCurrentValue();
}

G4String G4UImanager::GetCurrentStringValue(const char * aCommand, 
int parameterNumber, G4bool reGet)
{
  if(reGet || savedCommand == NULL)
  {
    savedParameters = GetCurrentValues( aCommand );
  }
  G4Tokenizer savedToken( savedParameters );
  G4String token;
  for(int i_thParameter=0;i_thParameter<parameterNumber;i_thParameter++)
  {
    token = savedToken();
    if( token.isNull() ) return G4String();
    if( token[(size_t)0] == '"' )
    {
      token.append(" ");
      token.append(savedToken("\""));
    }
  }
  return token;
}

G4String G4UImanager::GetCurrentStringValue(const char * aCommand, 
const char * aParameterName, G4bool reGet)
{
  if(reGet || savedCommand == NULL)
  {
    G4String parameterValues = GetCurrentValues( aCommand );
  }
  for(int i=0;i<savedCommand->GetParameterEntries();i++)
  {
    if( aParameterName ==
      savedCommand->GetParameter(i)->GetParameterName() )
      return GetCurrentStringValue(aCommand,i+1,false);
  }
  return G4String();
}

G4int G4UImanager::GetCurrentIntValue(const char * aCommand,
const char * aParameterName, G4bool reGet)
{
  G4String targetParameter = 
     GetCurrentStringValue( aCommand, aParameterName, reGet );
  G4int value;
  const char* t = targetParameter;
  G4std::istrstream is((char*)t);
  is >> value;
  return value;
}

G4int G4UImanager::GetCurrentIntValue(const char * aCommand,
int parameterNumber, G4bool reGet)
{
  G4String targetParameter = 
     GetCurrentStringValue( aCommand, parameterNumber, reGet );
  G4int value;
  const char* t = targetParameter;
  G4std::istrstream is((char*)t);
  is >> value;
  return value;
}

G4double G4UImanager::GetCurrentDoubleValue(const char * aCommand,
const char * aParameterName, G4bool reGet)
{
  G4String targetParameter = 
     GetCurrentStringValue( aCommand, aParameterName, reGet );
  G4double value;
  const char* t = targetParameter;
  G4std::istrstream is((char*)t);
  is >> value;
  return value;
}

G4double G4UImanager::GetCurrentDoubleValue(const char * aCommand,
int parameterNumber, G4bool reGet)
{
  G4String targetParameter = 
     GetCurrentStringValue( aCommand, parameterNumber, reGet );
  G4double value;
  const char* t = targetParameter;
  G4std::istrstream is((char*)t);
  is >> value;
  return value;
}

void G4UImanager::AddNewCommand(G4UIcommand * newCommand)
{
  treeTop->AddNewCommand( newCommand );
}

void G4UImanager::RemoveCommand(G4UIcommand * aCommand)
{
  treeTop->RemoveCommand( aCommand );
}

void G4UImanager::ExecuteMacroFile(const char * fileName)
{
  G4String fn = fileName;
  ExecuteMacroFile(fn);
}

void G4UImanager::ExecuteMacroFile(G4String fileName)
{
  G4UIsession* batchSession = new G4UIbatch(fileName,session);
  session = batchSession;
  G4UIsession* previousSession = session->SessionStart();
  delete session;
  session = previousSession;
}

void G4UImanager::LoopS(G4String valueList)
{
  G4Tokenizer parameterToken( valueList );
  G4String mf = parameterToken();
  G4String vn = parameterToken();
  G4String c1 = parameterToken();
  c1 += " ";
  c1 += parameterToken();
  c1 += " ";
  c1 += parameterToken();
  const char* t1 = c1;
  G4std::istrstream is((char*)t1);
  G4double d1;
  G4double d2;
  G4double d3;
  is >> d1 >> d2 >> d3;
  Loop(mf,vn,d1,d2,d3);
}
  
void G4UImanager::Loop(const char * macroFile,const char * variableName,
                   G4double initialValue,G4double finalValue,G4double stepSize)
{
  G4String mf = macroFile;
  G4String vn = variableName;
  Loop(mf,vn,initialValue,finalValue,stepSize);
}

void G4UImanager::Loop(G4String macroFile,G4String variableName,
                   G4double initialValue,G4double finalValue,G4double stepSize)
{
  G4String cd;
  for(G4double d=initialValue;d<=finalValue;d+=stepSize)
  {
    char st[20];
    G4std::ostrstream os(st,20);
    os << d << '\0';
    cd += st;
    cd += " ";
  }
  Foreach(macroFile,variableName,cd);
}

void G4UImanager::ForeachS(G4String valueList)
{
  G4Tokenizer parameterToken( valueList );
  G4String mf = parameterToken();
  G4String vn = parameterToken();
  G4String c1 = parameterToken();
  G4String ca;
  while(!((ca=parameterToken()).isNull()))
  {
    c1 += " ";
    c1 += ca;
  }
  Foreach(mf,vn,c1);
}

void G4UImanager::Foreach(const char * macroFile,const char * variableName,
                   const char * candidates)
{
  G4String mf = macroFile;
  G4String vn = variableName;
  G4String cd = candidates;
  Foreach(mf,vn,cd);
}

void G4UImanager::Foreach(G4String macroFile,G4String variableName,G4String candidates)
{
  G4Tokenizer parameterToken( candidates );
  G4String cd;
  while(!((cd=parameterToken()).isNull()))
  {
    SetAlias(variableName+" "+cd);
    ExecuteMacroFile(macroFile);
  }
}


int G4UImanager::ApplyCommand(const char * aCommand)
{
  G4String theCommand = aCommand;
  return ApplyCommand(theCommand);
}

G4String G4UImanager::SolveAlias(G4String aCmd)
{
  G4String aCommand = aCmd;
  int ia = aCommand.index("{");
  int iz = aCommand.index("#");
  while((ia != int(G4std::string::npos))&&((iz==int(G4std::string::npos))||(ia<iz)))
  {
    int ibx = -1;
    while(ibx<0)
    {
      int ib = aCommand.index("}");
      if( ib == int(G4std::string::npos) )
      {
        G4cerr << aCommand << G4endl;
        for(int iz=0;iz<ia;iz++) G4cerr << " ";
        G4cerr << "^" << G4endl;
        G4cerr << "Unmatched alias parenthis -- command ignored" << G4endl;
        return fAliasNotFound;
      }
      G4String ps = aCommand(ia+1,aCommand.length()-(ia+1));
      int ic = ps.index("{");
      int id = ps.index("}");
      if(ic!=int(G4std::string::npos) && ic < id)
      { ia+=ic+1; }
      else
      { ibx = ib; }
    }
    //--- Here ia represents the position of innermost "{"
    //--- and ibx represents corresponding "}"
    G4String subs;
    if(ia>0) subs = aCommand(0,ia);
    G4String alis = aCommand(ia+1,ibx-ia-1);
    G4String rems = aCommand(ibx+1,aCommand.length()-ibx);
    // G4cout << "<" << subs << "> <" << alis << "> <" << rems << ">" << G4endl;
    G4String* alVal = aliasList->FindAlias(alis);
    if(!alVal)
    {
      G4cerr << "Alias <" << alis << "> not found -- command ignored" << G4endl;
      G4String nullStr;
      return nullStr;
    }
    aCommand = subs+(*alVal)+rems;
    ia = aCommand.index("{");
  }
  return aCommand;
}

int G4UImanager::ApplyCommand(G4String aCommand)
{
  G4String aCmd = SolveAlias(aCommand);
  if(aCmd.isNull()) return fAliasNotFound;
  aCommand = aCmd;
  if(verboseLevel) G4cout << aCommand << G4endl;
  G4String commandString;
  G4String commandParameter;

  int i = aCommand.index(" ");
  if( i != int(G4std::string::npos) )
  {
    commandString = aCommand(0,i);
    commandParameter = aCommand(i+1,aCommand.length()-(i+1));
  }
  else
  {
    commandString = aCommand;
  }

  G4UIcommand * targetCommand = treeTop->FindPath( commandString );
  if( targetCommand == NULL )
  { return fCommandNotFound; }

  if(!(targetCommand->IsAvailable())) 
  { return fIllegalApplicationState; }
 
  if(saveHistory) historyFile << aCommand << G4endl; 
  histVec.push_back(aCommand);
  return targetCommand->DoIt( commandParameter );
}

void G4UImanager::StoreHistory(G4String fileName)
{ StoreHistory(true,fileName); }

void G4UImanager::StoreHistory(G4bool historySwitch,G4String fileName)
{
  if(historySwitch)
  {
    if(saveHistory)
    { historyFile.close(); }
    const char* theFileName = fileName;
    historyFile.open((char*)theFileName);
    saveHistory = true;
  }
  else
  {
    historyFile.close();
    saveHistory = false;
  }
  saveHistory = historySwitch;
}

void G4UImanager::PauseSession(G4String msg)
{
  if(session) session->PauseSessionStart(msg); 
}

void G4UImanager::ListCommands(G4String direct)
{
  G4UIcommandTree* comTree = FindDirectory(direct);
  if(comTree)
  { comTree->List(); }
  else
  { G4cout << direct << " is not found." << G4endl; }
}

G4UIcommandTree* G4UImanager::FindDirectory(const char* dirName)
{
  G4String aDirName = dirName;
  G4String targetDir = aDirName.strip(G4String::both);
  if( targetDir( targetDir.length()-1 ) != '/' )
  { targetDir += "/"; }
  G4UIcommandTree* comTree = treeTop;
  if( targetDir == "/" )
  { return comTree; }
  int idx = 1;
  while( idx < int(targetDir.length())-1 )
  {
    int i = targetDir.index("/",idx);
    comTree = comTree->GetTree(targetDir(0,i+1));
    if( comTree == NULL )
    { return NULL; }
    idx = i+1;
  }
  return comTree;
}

G4bool G4UImanager::Notify(G4ApplicationState requestedState)
{
  //G4cout << G4StateManager::GetStateManager()->GetStateString(requestedState) << " <--- " << G4StateManager::GetStateManager()->GetStateString(G4StateManager::GetStateManager()->GetPreviousState()) << G4endl;
  if(pauseAtBeginOfEvent)
  {
    if(requestedState==EventProc &&
       G4StateManager::GetStateManager()->GetPreviousState()==GeomClosed)
    { PauseSession("BeginOfEvent"); }
  }
  if(pauseAtEndOfEvent)
  {
    if(requestedState==GeomClosed &&
       G4StateManager::GetStateManager()->GetPreviousState()==EventProc)
    { PauseSession("EndOfEvent"); }
  }
  return true;
}

void G4UImanager::Interact()
{
  Interact(G4String("G4> "));
}

void G4UImanager::Interact(const char * pC)
{
  G4String cc = pC;
  Interact(cc);
}

void G4UImanager::Interact(G4String pC)
{
  G4cerr << "G4UImanager::Interact() is out of date and is not used anymore." << G4endl;
  G4cerr << "This method will be removed shortly!!!" << G4endl;
  G4cerr << "In case of main() use" << G4endl;
  G4cerr << "    G4UIsession * session = new G4UIterminal;" << G4endl;
  G4cerr << "    session->SessionStart();" << G4endl;
  G4cerr << "In other cases use" << G4endl;
  G4cerr << "    G4StateManager::GetStateManager()->Pause();" << G4endl;
}



void G4UImanager::SetCoutDestination(G4UIsession *const value)
{
    G4coutbuf.SetDestination(value);
    G4cerrbuf.SetDestination(value);
}

void G4UImanager::SetAlias(const char * aliasLine)
{
  G4String aL = aliasLine;
  SetAlias(aL);
}

void G4UImanager::SetAlias(G4String aliasLine)
{
  int i = aliasLine.index(" ");
  G4String aliasName = aliasLine(0,i);
  G4String aliasValue = aliasLine(i+1,aliasLine.length()-(i+1));
  if(aliasValue(0)=='"')
  { 
    G4String strippedValue;
    if(aliasValue(aliasValue.length()-1)=='"')
    { strippedValue = aliasValue(1,aliasValue.length()-2); }
    else
    { strippedValue = aliasValue(1,aliasValue.length()-1); }
    aliasValue = strippedValue;
  }

  aliasList->ChangeAlias(aliasName,aliasValue);
}

void G4UImanager::RemoveAlias(const char * aliasName)
{
  G4String aL = aliasName;
  RemoveAlias(aL);
}

void G4UImanager::RemoveAlias(G4String aliasName)
{
  G4String targetAlias = aliasName.strip(G4String::both);
  aliasList->RemoveAlias(targetAlias);
}

void G4UImanager::ListAlias()
{
  aliasList->List();
}

void G4UImanager::CreateHTML(const char* dir)
{ 
  G4UIcommandTree* tr = FindDirectory(dir);
  if(tr!=0)
  { tr->CreateHTML(); }
  else
  { G4cerr << "Directory <" << dir << "> is not found." << G4endl; }
}

