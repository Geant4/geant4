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
// 12/06/2002 G4UIGainServer H. MInamimoto and H. Yoshida created
// $Id: G4UIGainServer.cc 66892 2013-01-17 10:57:59Z gunter $
//
#ifndef WIN32

#include "G4UIGainServer.hh"
#include <netdb.h>

#include <sstream>
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"


//////////////////////////////////////////////
G4UIGainServer::G4UIGainServer()
//////////////////////////////////////////////
{
    TVersion ="T1.0a"; JVersion="J1.0a";
    prefix = "/";

    port = DEFAULT_PORT;
    while(SetUPServer() == false){
        G4cout<<"can't get the port no. "<<port<<" Now, try to get the next port "<<port+1<<G4endl;
        port++;
    }


    UI= G4UImanager::GetUIpointer();
    UI-> SetSession(this);
    UI-> SetCoutDestination(this);

    G4StateManager* statM = G4StateManager::GetStateManager();
    promptCharacter = statM->GetStateString(statM->GetCurrentState());
    uiMode = terminal_mode;

    iExit= FALSE;
    iCont= FALSE;

    G4UIcommandTree* tree = UI->GetTree();
    GetNewTreeStructure(tree,0);
    GetNewTreeValues(tree,0);
    previousTreeCommands = newTreeCommands;
    previousTreeParams = newTreeParams;
    previousTreePCP = newTreePCP;

}

/////////////////////////////
G4UIGainServer::~G4UIGainServer() 
/////////////////////////////
{ 

    if(G4UImanager::GetUIpointer()) {
      UI-> SetSession(NULL);
      UI-> SetCoutDestination(NULL);
    }

    if(G4UImanager::GetUIpointer()!=0){
        UI->SetSession(NULL);
        UI->SetCoutDestination(NULL);
    }
}


/////////////////////////////////////////
G4UIsession* G4UIGainServer::SessionStart()
/////////////////////////////////////////
{
    G4String newCommand;

    G4StateManager* statM = G4StateManager::GetStateManager();
    promptCharacter = statM->GetStateString(statM->GetCurrentState());
    
    iExit= TRUE;

    WaitingConnection();
    while(iExit){
        newCommand= GetCommand();
        ExecuteCommand(newCommand);
    }
    return NULL;
}

//////////////////////////////////////////////////
void G4UIGainServer::PauseSessionStart(const G4String& msg)
//////////////////////////////////////////////////
{
    promptCharacter = msg;
    G4cout<<"@@PROMPT \""<<promptCharacter<<"\""<<G4endl;

    iCont= TRUE;

    G4String newCommand= GetCommand();
    while(iCont){
      ExecuteCommand(newCommand);
      newCommand= GetCommand();
      strcpy(buf,"nowIdle");
      write(socketD[1],buf,strlen(buf));
    }
}

////////////////////////////////////////////////////
void G4UIGainServer::ExecuteCommand(const G4String& aCommand)
////////////////////////////////////////////////////
{
    if(aCommand.length()<2) return;

    G4UIcommandTree* tree = UI->GetTree();
    if(aCommand.length()<2) return;
    G4int returnVal = UI->ApplyCommand(aCommand);
    G4int paramIndex = returnVal % 100;
    // 0 - 98 : paramIndex-th parameter is invalid
    // 99     : convination of parameters is invalid
    G4int commandStatus = returnVal - paramIndex;

    UpdateState();

    if(uiMode != terminal_mode){
        switch(commandStatus) {
        case fCommandSucceeded:
            GetNewTreeStructure(tree,0);
            GetNewTreeValues(tree,0);
            if(CommandUpdated()){
                NotifyCommandUpdate();
            } else{
                UpdateParamVal();
            }
            previousTreeCommands = newTreeCommands;
            previousTreeParams = newTreeParams;
            previousTreePCP = newTreePCP;
            break;
        case fCommandNotFound:
            G4cerr << "@@ErrResult \" <" << UI->SolveAlias(aCommand) << "> not found.\"" << G4endl;
            break;
        case fIllegalApplicationState:
            G4cerr << "@@ErrResult \"illegal application state -- command refused.\"" << G4endl;
            break;
        case fParameterOutOfRange:
            G4cout << "@@ErrResult \"Parameter Out of Range.\"" << G4endl;
            break;
        case fParameterUnreadable:
            G4cout << "@@ErrResult \"Parameter is wrong type and/or is not omittable.\""<<G4endl;
            break;
        case fParameterOutOfCandidates:
            G4cerr << "@@ErrResult \"Parameter is out of candidate.\"" << G4endl;
            break;
        case fAliasNotFound:
        default:
            G4cerr << "command refused (" << commandStatus << ")" << G4endl;
        }
    }
}

///////////////////////////////////
G4String G4UIGainServer::GetCommand()
///////////////////////////////////
{
    G4String newCommand;
    G4String nullString;

  while( 1 )
  {
    G4UIcommandTree* tree = UI->GetTree();
    if ( uiMode != terminal_mode ){
      G4cout << "@@PROMPT \"" << promptCharacter << "\"" << G4endl;
    }
    if ( uiMode != java_mode ){
      G4cout << promptCharacter << "> " << G4endl;
    }else{
      G4cout << "@@Ready" << G4endl;
    }


    /////////////////////////////
    /////////////////////////////
    read(socketD[1],buf,1024);
    newCommand=buf;
    //DEBUG cout<<"->"<<newCommand<<"<-"<<newCommand.length()<<G4endl;
    //newCommand.readLine( G4cin, FALSE );
    /////////////////////////////
    /////////////////////////////



    if (!G4cin.good()) { G4cin.clear(); newCommand = nullString; iExit=false;break;}

    newCommand = newCommand.strip(G4String::leading);
    if( newCommand.length() < 1) { break; }

    while( newCommand(newCommand.length()-1) == '_' )
    {
      G4String newLine;
      newCommand.remove(newCommand.length()-1);
      newLine.readLine( G4cin );
      if (!G4cin.good()) { G4cin.clear(); newCommand = nullString; iExit=false;break;}
      newCommand.append(newLine);
    }

    G4String nC = newCommand.strip(G4String::leading);
    if( nC.length() < 1) { break; }

    // -------------------- nC.toUpper();
    if( nC == "@@GainmodeJAVA" ) {
      uiMode = java_mode;
      G4cout << G4endl << "@@Version " << JVersion << G4endl;
      SendCommandProperties(tree);
      NotifyStateChange();
    }
    else if( nC == "@@GainmodeTcl" ) {
      uiMode = tcl_mode;
      G4cout << G4endl << "@@Version " << TVersion << G4endl;
      SendCommandProperties(tree);
      NotifyStateChange();
    }
    else if( nC(0) == '#' )
      { G4cout << nC << G4endl; }

    else if( nC == "ls"  || nC(0,3) == "ls " )
    { ListDirectory( nC ); }
    else if( nC == "pwd" )
    { G4cout << "Current Working Directory : " << prefix << G4endl; }
    else if( nC(0,2) == "cd"  || nC(0,3) == "cd " )
    { ChangeDirectory( nC ); }
    else if(  nC == "help" || nC(0,5) == "help ")
    { TerminalHelp( nC ); }
    else if( nC(0) == '?' )
    { ShowCurrent( nC ); }
    else if( nC(0,4) == "hist"   || nC == "history")
    {
      G4int nh = UI->GetNumberOfHistory();
      for(int i=0;i<nh;i++)
      { G4cout << i << ": " << UI->GetPreviousCommand(i) << G4endl; }
    }
    else if( nC(0) == '!' )
    {
      G4String ss = nC(1,nC.length()-1);
      G4int vl;
      const char* tt = ss;
      std::istringstream is((char*)tt);
      is >> vl;
      G4int nh = UI->GetNumberOfHistory();
      if(vl>=0 && vl<nh)
      {
        newCommand = UI->GetPreviousCommand(vl);
        G4cout << newCommand << G4endl;
        break;
      }
      else
      { G4cerr << "history " << vl << " is not found." << G4endl; }
    }
    else if( nC(0,4) == "exit" )
    {
      if( iCont )
      {
        if ( uiMode == terminal_mode){
          G4cerr << "You are now processing RUN." << G4endl;
          G4cerr << "Please abrot it using \"/run/abort\" command first" << G4endl;
          G4cerr << " and use \"continue\" command until the application" << G4endl;
          G4cerr << " becomes to Idle." << G4endl;
        }else{
          G4cout << "@@ErrResult \"You are now processing RUN.\"" << G4endl;
        }
      }
      else
      {
        close(socketD[1]);
        close(socketD[2]);
        iExit = false;
        newCommand = nullString;
        break;
      }
    }
    else if(  nC == "cont" || nC == "continue" )
    {
      iCont = false;
      newCommand = nullString;
      break;
    }
    else
    { break; }
  }
  return GetFullPath(newCommand);
}

//////////////////////////////////////////////////////
G4int G4UIGainServer::ReceiveG4cout(const G4String& coutString)
//////////////////////////////////////////////////////
{
    if(socketD[1]>0){
        write(socketD[1],coutString,coutString.length());
    }
    return 0;

  //std::cout << coutString << std::flush;
  //return 0;
}

//////////////////////////////////////////////////////
G4int G4UIGainServer::ReceiveG4cerr(const G4String& cerrString)
//////////////////////////////////////////////////////
{
    if(socketD[2]>0){
        write(socketD[2],cerrString,cerrString.length());
    }
    return 0;

  //std::cerr << cerrString << std::flush;
  //return 0;
}

/////////////////////////////////////////////////
G4bool G4UIGainServer::GetHelpChoice(G4int& aInt)
/////////////////////////////////////////////////
{
    G4cin >> aInt;
    if(!G4cin.good()){
        G4cin.clear();
        G4cin.ignore(30,'\n');
        return FALSE;
    }
    return TRUE;
}

/////////////////////////////
void G4UIGainServer::ExitHelp() const
/////////////////////////////
{
    char temp[100];
    G4cin.getline(temp, 100);
}

/////////////////////////////
bool G4UIGainServer::SetUPServer(){
/////////////////////////////

    socketD[0] = socket(AF_INET,SOCK_STREAM,0);

    if(socketD[0]<0){
        perror("server:socket");
        return (false);
        //exit(1);
    }

    memset( (char *)&saddr,'\0',sizeof(saddr)) ;

    saddr.sin_family = AF_INET;
    saddr.sin_addr.s_addr = INADDR_ANY;
    saddr.sin_port = htons(port);
    unlink(SOCK_NAME);    

    if(bind(socketD[0] , (struct sockaddr *)&saddr , sizeof(saddr))<0){
        perror("bind");
        return (false);
        //exit(1);
    }
    else{ G4cout<<"G4GainServer waiting at "<<port<<G4endl; }

    if(listen(socketD[0],1)<0){
        perror("listen");
        return (false);
        //exit(1);
    }

    return (true);
}

////////////////////////////////////////
void G4UIGainServer::WaitingConnection(){
////////////////////////////////////////
    len = sizeof(caddr);

    for(int i=1;i<=2;i++){
#if defined __APPLE__ && (__GNUC__<4)
        if((socketD[i] = accept(socketD[0], (struct sockaddr *)&caddr,(int *)&len))<0){
#else
        if((socketD[i] = accept(socketD[0], (struct sockaddr *)&caddr,(socklen_t *)&len))<0){
#endif
            G4cerr<<"accept:"<<i<<G4endl;
            //exit(1);
            G4Exception("G4UIGainServer::WaitingConnection()",
                        "UI0004",
                        FatalException,
                        "Invalid Socket. Cannot establish connection");
        }
    }
    close(socketD[0]);
}

///////////////////////////////////////////////////
G4String G4UIGainServer::GetFullPath(G4String aNewCommand){
///////////////////////////////////////////////////
  G4String newCommand = aNewCommand.strip(G4String::both);
  G4String tmpString;
  if( newCommand(0) == '/' ) 
  { tmpString = newCommand; }
  else if( newCommand(0,3) == "../" )
  {
    G4String tmpPrefix = prefix;
    /*G4int*/ unsigned i_direc = 0;
    while( i_direc < newCommand.length() )
    { 
      if( newCommand(i_direc,3) == "../" )
      {
        i_direc += 3;
        prefix = ModifyPrefix( G4String("../") );
      }
      else
      { break; }
    }
    tmpString = prefix;
    tmpString.append( newCommand( i_direc, newCommand.length()-i_direc ) );
    prefix = tmpPrefix;
  }
  else
  {
    tmpString = prefix;
    tmpString.append( newCommand );
  }
  return tmpString;
}

////////////////////////////////
void G4UIGainServer::SessionTerminate(){
////////////////////////////////
    G4cout<<"***** Terminal session end *****"<<G4endl;
}


//////////////////////////////////////////////
void G4UIGainServer::ShowCurrent(G4String newCommand){
//////////////////////////////////////////////
  G4String theCommand = GetFullPath(newCommand(1,newCommand.length()-1));
  G4String curV = UI->GetCurrentValues(theCommand);
  if( ! (curV.isNull()||curV(0)=='\0' ) ) {
    if (uiMode == terminal_mode){
      G4cout << "Current value(s) of the parameter(s) : " << curV << G4endl;
    }else{
      G4cout << "@@CurrentValue " << curV << G4endl;
    }
  } else if (uiMode == terminal_mode){
      G4cout << "Current value is not available." << G4endl;
    } else {
      G4cout << "@@ErrResult \"Current value is not available.\"" << G4endl;
    }
}

//////////////////////////////////////////////////
void G4UIGainServer::ChangeDirectory(G4String newCommand){
//////////////////////////////////////////////////
  G4String savedPrefix = prefix;
  if( newCommand.length() <= 3 )
  { prefix = "/"; }
  else
  { 
    G4String aNewPrefix = newCommand(3,newCommand.length()-3);
    G4String newPrefix = aNewPrefix.strip(G4String::both);
    if( newPrefix(0) == '/' )
    { prefix = newPrefix; }
    else if( newPrefix(0) != '.' )
    { 
      prefix += newPrefix;
    }
    else
    { prefix = ModifyPrefix( newPrefix ); }
  }
  if( prefix( prefix.length() - 1 ) != '/' )
  { prefix += "/"; }
  if( FindDirPath( prefix ) == NULL )
  {
    G4cout << "Directory <" << prefix << "> is not found." << G4endl;
    prefix = savedPrefix;
  }
}
////////////////////////////////////////////////
void G4UIGainServer::ListDirectory(G4String newCommand){
////////////////////////////////////////////////
  G4String targetDir('\0');
  if( newCommand.length() <= 3 )
  { targetDir = prefix; }
  else
  {
    G4String newPrefix = newCommand(3,newCommand.length()-3);
    newPrefix.strip(G4String::both);
    if( newPrefix(0) == '/' )
    { targetDir = newPrefix; }
    else if( newPrefix(0) != '.' )
    {
      targetDir = prefix;
      targetDir += newPrefix;
    }
    else
    { targetDir = ModifyPrefix( newPrefix ); }
  }
  if( targetDir( targetDir.length() - 1 ) != '/' )
  { targetDir += "/"; }
  G4UIcommandTree * commandTree = FindDirPath( targetDir );
  if( commandTree == NULL )
  { G4cout << "Directory <" << targetDir << "> is not found." << G4endl; }
  else
  { commandTree->ListCurrent(); }
}

///////////////////////////////////////////////
void G4UIGainServer::TerminalHelp(G4String newCommand){
///////////////////////////////////////////////
    G4UIcommandTree* treeTop = UI->GetTree();
    str_size i = newCommand.index(" ");
    
    if(i!=std::string::npos){
        G4String newValue = newCommand(i+1,newCommand.length()-(i+1));
        newValue.strip(G4String::both);
        if(newValue(0)!='/'){
            newValue.prepend(prefix);
        }
        G4UIcommand* theCommand = treeTop->FindPath(newValue);
        if(theCommand !=NULL){
            theCommand->List();
            return;
        }
        else{
            G4cout<<"Command<" << newValue << "is not found."<<G4endl;
            return;
        }
    }

    G4UIcommandTree* floor[10];
    floor[0] = treeTop;
    int iFloor = 0;
    unsigned prefixIndex = 1;
    while(prefixIndex<prefix.length()-1){
        int ii = prefix.index("/",prefixIndex);
        floor[iFloor+1]=
          floor[iFloor]->GetTree(G4String(prefix(0,ii+1)));
        prefixIndex = ii+1;
        iFloor++;
    }
    floor[iFloor]->ListCurrentWithNum();
    while(1){
        int j;
        G4cout<<G4endl <<"Type the number (0:end, -n:n level back) :"<<std::flush;
        G4cin >> j;
        if(!G4cin.good()){
            G4cin.clear();
            G4cin.ignore(30,'\n');
            G4cout<<G4endl <<"Not a number,once more"<<G4endl; continue;
        }
        else if(j<0){
            iFloor += j;
            if(iFloor <0) iFloor =0;
            floor[iFloor]->ListCurrentWithNum(); continue;
        }
        else if(j==0){break;}
        else if(j>0){
            int n_tree = floor[iFloor]->GetTreeEntry();
            if(j>n_tree){
                if(j<=n_tree+floor[iFloor]->GetCommandEntry()){
                    floor[iFloor]->GetCommand(j-n_tree)->List();
                }
            }
            else{
                floor[iFloor+1] = floor[iFloor]->GetTree(j);
                iFloor++;
                floor[iFloor]->ListCurrentWithNum();
            }
        }
    }
    G4cout<<"Exit from Help."<<G4endl <<G4endl;
    G4cout<<G4endl;
    char temp[100];
    G4cin.getline(temp,100);
}


///////////////////////////////////////////////////
G4String G4UIGainServer::ModifyPrefix(G4String newCommand){
///////////////////////////////////////////////////
    G4String newPrefix = prefix;
    while(1){
        if(newCommand(0,2) ==".."){
            if(newPrefix !="/"){
                G4String tmpString = newPrefix(0,newPrefix.length()-1);
                newPrefix = newPrefix(0,tmpString.last('/')+1);
            }
        }
        else{
            newPrefix += newCommand;
            break;
        }
        if(newCommand == ".." || newCommand == "../"){
            break;
        }
        newCommand=newCommand(3,newCommand.length()-3);
    }
    return newPrefix;
}

//////////////////////////////////////////////////////////
G4UIcommandTree* G4UIGainServer::FindDirPath(G4String newCommand){
//////////////////////////////////////////////////////////
  G4UIcommandTree * comTree = UI->GetTree();
  /*int*/ unsigned idx = 1; 
  while( idx < newCommand.length()-1 )
  { 
    int i = newCommand.index("/",idx);
    comTree = comTree->GetTree(G4String(newCommand(0,i+1)));
    if( comTree == NULL )
    { return NULL; }
    idx = i+1;
  }
  return comTree;
}

//// ----- for JAVA Gain

//////////////////////////////////////////////////////////
void G4UIGainServer::SendCommandProperties(G4UIcommandTree* tree){
//////////////////////////////////////////////////////////
  if( tree == NULL ) {
    G4cerr << "GetTree() returnes null." << G4endl;
    return;
  }
  if (uiMode == java_mode){
    G4cout << "@@JTreeBegin" << G4endl;
    CodeGenJavaTree(tree, 0);
    G4cout << "@@JTreeEnd" << G4endl;
    CodeGenJavaParams(tree, 0);
  }else{}
}

////////////////////////////////////////////////////////////
void G4UIGainServer::SendParameterProperties(G4UIcommandTree* tree){
////////////////////////////////////////////////////////////
  if( tree == NULL ) {
    G4cerr << "GetTree() returnes null." << G4endl;
    return;
  }
  if (uiMode == java_mode){
    CodeGenJavaParams(tree, 0);
  }else{ }
}

//////////////////////////////////////////////////////////////
void G4UIGainServer::CodeGenJavaTree(G4UIcommandTree* tree,int level){
//////////////////////////////////////////////////////////////
  int treeEntry, commandEntry;
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  if(level!=0) {
    for(int i=0; i<commandEntry; i++){
      G4cout << tree->GetCommand(i+1)->GetCommandPath() << G4endl;
    }
  }
  if(treeEntry == 0) return; //end recursion

  for(int j=0; j<treeEntry; j++){
    CodeGenJavaTree(tree->GetTree(j+1), level+1);
  }
}

////////////////////////////////////////////////////////////////
void G4UIGainServer::CodeGenJavaParams(G4UIcommandTree* tree,int level){
////////////////////////////////////////////////////////////////
    int treeEntry,commandEntry,i;
    G4UIcommandTree* treeLink;

    treeEntry = tree->GetTreeEntry();
    commandEntry = tree->GetCommandEntry();

    for(i=0;i<commandEntry; i++){
        SendAParamProperty(tree->GetCommand(i+1));
    }
    if(treeEntry ==0) return;

    for(i=0;i<treeEntry; i++){
        treeLink = tree->GetTree(i+1);
        G4cout<<"@@JDirGuieBegin"<<G4endl;
        G4cout<<treeLink->GetPathName()<<G4endl <<treeLink->GetTitle()<<G4endl;
        G4cout<<"@@JDirGuideEnd"<<G4endl;
        CodeGenJavaParams(treeLink,level+1);
    }
}

///////////////////////////////////////////////////
void G4UIGainServer::SendAParamProperty(G4UIcommand* Comp){
///////////////////////////////////////////////////
  int guidanceEntry, parameterEntry;
  G4String title, title2;
  G4UIparameter * prp;
  char c[2];
  guidanceEntry = Comp->GetGuidanceEntries();
  parameterEntry = Comp->GetParameterEntries();
  G4cout << "@@JParamBegin" << G4endl;
  G4cout << Comp->GetCommandPath() << G4endl;
  G4cout << guidanceEntry << G4endl;
  for (int j=0; j<guidanceEntry; j++){
    title = Comp->GetGuidanceLine(j);
    title2 = "";
    if (title != ""){
      for(int i=0; i< (int)title.length(); i++){
        c[0]=title(i);
        c[1]= '\0';
        if ( c[0] == '\n' || c[0] == '\r') {
          c[0]= ' ';
        }
        title2.append(c);
      }
    }
    G4cout << title2 << G4endl;
  }
  G4cout << Comp->GetRange() << G4endl;
  G4cout << parameterEntry << G4endl;
  for( int par=0; par<parameterEntry; par++) {
    prp = (G4UIparameter *)Comp->GetParameter(par);
    G4cout << prp->GetParameterName() << G4endl;
    G4cout << prp->GetParameterGuidance() << G4endl;
    G4cout << prp->GetParameterType() << G4endl;
    G4cout << prp->IsOmittable() << G4endl;
    G4cout << prp->GetDefaultValue() << G4endl;
    G4cout << prp->GetParameterRange() << G4endl;
    G4cout << prp->GetParameterCandidates() << G4endl;
  }
  G4cout << "@@JParamEnd" << G4endl;
}

//////////////////////////////////////////////////////////////
void G4UIGainServer::SendDisableList(G4UIcommandTree* tree,int level){
//////////////////////////////////////////////////////////////
  int treeEntry, commandEntry;
  G4UIcommand * Comp;
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  for(int com=0; com<commandEntry; com++) {
    Comp = tree->GetCommand(com+1);
    if( Comp->IsAvailable()==false ) {
       G4cout << Comp->GetCommandPath()<<G4endl;
    }
  }
  if( treeEntry == 0 ) return;     // end recursion

  for( int i=0; i<treeEntry; i++) {
    SendDisableList(tree->GetTree(i+1), level+1);
    // be sure the function name is the same
  }
}



//####### update check routines ####################################

///////////////////////////////
void G4UIGainServer::UpdateState(void)
///////////////////////////////
{
   static G4ThreadLocal G4ApplicationState *previousState_G4MT_TLS_ = 0 ; if (!previousState_G4MT_TLS_) {previousState_G4MT_TLS_ = new  G4ApplicationState  ; *previousState_G4MT_TLS_= G4State_PreInit ; }  G4ApplicationState &previousState = *previousState_G4MT_TLS_;
   G4ApplicationState  newState;
   G4StateManager *statM = G4StateManager::GetStateManager();
   newState = statM->GetCurrentState();
   if( newState != previousState ) 
   {
      NotifyStateChange();
      previousState = newState; 
   }
}

/////////////////////////////////////
void G4UIGainServer::NotifyStateChange(void)
/////////////////////////////////////
{
   G4String stateString;
   G4StateManager * statM = G4StateManager::GetStateManager();
   G4UIcommandTree * tree = UI->GetTree();
   stateString = statM->GetStateString(statM->GetCurrentState());
   if ( uiMode != terminal_mode ){
     G4cout << "@@State \"" << stateString << "\"" << G4endl;
     G4cout << "@@DisableListBegin"<<G4endl;
     SendDisableList(tree, 0);
     G4cout << "@@DisableListEnd" <<G4endl;
   }
}

///////////////////////////////////////
void G4UIGainServer::NotifyCommandUpdate(void)
///////////////////////////////////////
{
  G4UIcommandTree * tree = UI->GetTree();
  SendCommandProperties(tree);
}

/////////////////////////////////////////////////////
void G4UIGainServer::NotifyParameterUpdate(G4UIcommand* com)
/////////////////////////////////////////////////////
{
    SendAParamProperty(com);
}

//////////////////////////////////
int G4UIGainServer::CommandUpdated(void){
//////////////////////////////////
  int added=0, deleted=0;
  int pEntry= previousTreeCommands.size();
  int nEntry= newTreeCommands.size();
  int i,j;
  for( i=0; i<pEntry; i++) {      // check deleted command(s)
      for( j=0; j<nEntry; j++) {
         if( previousTreeCommands[i] == newTreeCommands[j]) break;
      }
      if( j==nEntry ) { 
         deleted = 1;
         //G4cout <<"deleted: "<< previousTreeCommands(i) << G4endl;
      }
  }
  for( i=0; i<nEntry; i++) {      // check added command(s)
      for( j=0; j<pEntry; j++) {
         if( newTreeCommands[i] == previousTreeCommands[j]) break;
      }
      if( j==pEntry ) { 
         added = 1;
      //   G4cout <<"added: "<< newTreeCommands(i) << G4endl;
      }
  }
  if( added    && deleted==0 ) {G4cout<<"c added"<<G4endl;return added;}
  if( added==0 && deleted ) {G4cout<<"c deleted"<<G4endl;return deleted;}
  if( added    && deleted ) {G4cout<<"c add/deleted"<<G4endl;return addedAndDeleted;}
  return notChanged;
}

//////////////////////////////////////////////////////////////////////
void G4UIGainServer::GetNewTreeStructure(G4UIcommandTree * tree, int level) { 
//////////////////////////////////////////////////////////////////////
  G4String commandPath;
  G4String title; 
  G4String pathName; //tree name
  G4UIcommandTree * t;
  int treeEntry    = tree->GetTreeEntry();
  int commandEntry = tree->GetCommandEntry();

  if( level==0 ) { newTreeCommands.clear();}
  for(int com=0; com<commandEntry; com++){
      commandPath = tree->GetCommand(com+1)->GetCommandPath();
      title = tree->GetCommand(com+1)->GetTitle();
      newTreeCommands.push_back( commandPath + " " + title );
  }

  if(treeEntry == 0) return; //end recursion

  for(int i=0; i< treeEntry; i++){
    t = tree->GetTree(i+1);
    pathName =  t->GetPathName();   
    title = t->GetTitle();
    newTreeCommands.push_back( pathName + " " + title );
    GetNewTreeStructure(t, level+1);
  }
}

////////////////////////////////////
void G4UIGainServer::UpdateParamVal(void) {
////////////////////////////////////
  // call NotifyParameterUpdate() if the value of each
  //  command/parameter is updated.
  //  assuming the command structure is not changed.
  int pEntry= previousTreeParams.size();
  int nEntry= newTreeParams.size();
  int i;
  G4UIcommand* Comp;
  if (pEntry != nEntry) return; 
  for( i=0; i<nEntry; i++) {
    if( previousTreeParams[i] != newTreeParams[i]){
       Comp = newTreePCP[i];
       G4cout << Comp->GetCommandPath()
            << " command is updated." <<G4endl; 
       NotifyParameterUpdate(Comp);
    }
  }
}

//////////////////////////////////////////////////////////////////
void G4UIGainServer::GetNewTreeValues( G4UIcommandTree * tree, int level){ // recursive
//////////////////////////////////////////////////////////////////
   G4String commandPath;
   G4String pathName; //tree name
   G4UIcommandTree * t;
   int parameterEntry;
   int treeEntry    = tree->GetTreeEntry();
   int commandEntry = tree->GetCommandEntry();
   G4UIcommand * Comp;
   G4UIparameter * prp; 
   G4String param, str(" ");

   if( level==0 ) { newTreeParams.clear(); }
   for(int com=0; com<commandEntry; com++) {
      Comp = tree->GetCommand(com+1);
      commandPath    = Comp->GetCommandPath();
      parameterEntry = Comp->GetParameterEntries();
      param = commandPath +" ";
      for( int par=0; par< parameterEntry; par++) {
         prp = (G4UIparameter *)Comp->GetParameter(par);
         param += prp->GetParameterName() +" ";
         str(0) = prp->GetParameterType();
         param += str + " ";
         param += prp->GetDefaultValue()  +" ";
         param += prp->GetParameterRange() +" ";
         param += prp->GetParameterCandidates();
      }
     newTreeParams.push_back( param + "\n"); 
     newTreePCP.push_back( Comp ); 
   }
   if( treeEntry == 0 )  return;     // end recursion
   for( int i=0; i< treeEntry; i++) {
      t = tree->GetTree(i+1);
      GetNewTreeValues(t, level+1);
   }
}


#endif

        









        
