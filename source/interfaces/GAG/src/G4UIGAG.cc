// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIGAG.cc,v 1.1 1999-01-07 16:09:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4UIGAG.cc
// 18.Feb.98 M.Nagamatu and T.Kodama created G4UIGAG from G4UIterminal

#include "G4UIGAG.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

#ifdef GNU_GCC
  template class RWTValVector<G4UIcommand *>;
  template class RWTValVector<RWCString>;
#endif

G4UIGAG::G4UIGAG(): TVersion("T1.0a"), JVersion("J1.0a")
{
  //G4cout << "G4UIGAG: Apr15,98." << endl;
  prefix = "/";
  UI = G4UImanager::GetUIpointer();
  UI->SetSession(this);
  UI->SetCoutDestination(this); 
  G4StateManager * statM = G4StateManager::GetStateManager();
  promptCharacter = statM->GetStateString(statM->GetCurrentState());
  uiMode = terminal_mode; // GAG
  iExit = false;
  iCont = false;
  // -- Initialize Notify routines begin
  G4UIcommandTree * tree = UI->GetTree();
  GetNewTreeStructure(tree,0);
  GetNewTreeValues(tree,0);
  previousTreeCommands = newTreeCommands;
  previousTreeParams = newTreeParams;
  previousTreePCP = newTreePCP;
  // -- end
}

G4UIGAG::~G4UIGAG()
{
  if( G4UImanager::GetUIpointer() != 0)
  {
     UI->SetSession(NULL);
     UI->SetCoutDestination(NULL);
     //     G4cout << "GAG session deleted" << endl;
  }
}                                                 

G4UIsession * G4UIGAG::SessionStart()
{
  iExit = true;
  G4StateManager * statM = G4StateManager::GetStateManager();
  promptCharacter = statM->GetStateString(statM->GetCurrentState());
  G4String newCommand = GetCommand();
  while( iExit )
  {
    ExecuteCommand(newCommand);
    promptCharacter = statM->GetStateString(statM->GetCurrentState());
    newCommand = GetCommand();
  }
  return NULL;
}

void G4UIGAG::PauseSessionStart(G4String msg)
{
  promptCharacter = msg;
  G4cout << "@@PROMPT \"" << promptCharacter << "\"" << endl;
  iCont = true;
  G4String newCommand = GetCommand();
  while( iCont )
  {
    ExecuteCommand(newCommand);
    newCommand = GetCommand();
  }
}

void G4UIGAG::ExecuteCommand(G4String aCommand)
{
  G4UIcommandTree * tree = UI->GetTree();
  if(aCommand.length()<2) return;
  int commandStatus = UI->ApplyCommand(aCommand);
  UpdateState();
  if ( uiMode == terminal_mode){
    switch(commandStatus) {
    case fCommandSucceeded:
      break;
    case fCommandNotFound:
      G4cerr << "command not found" << endl;
      break;
    case fIllegalApplicationState:
      G4cerr << "illegal application state -- command refused" << endl;
      break;
    case fParameterOutOfRange:
    case fParameterUnreadable:
    case fParameterOutOfCandidates:
    default:
      G4cerr << "command refused (" << commandStatus << ")" << endl;
    }
  }else{
    switch(commandStatus) {
    case fCommandSucceeded:
      {
        GetNewTreeStructure(tree,0);
        GetNewTreeValues(tree,0);
        if (CommandUpdated()) {
           NotifyCommandUpdate();
        } else {
           UpdateParamVal();  // if param is updated, call notifyPara...
        } 
        previousTreeCommands = newTreeCommands;
        previousTreeParams = newTreeParams;
        previousTreePCP = newTreePCP;
      }
      break;
    case fCommandNotFound:
      G4cout << "@@ErrResult \"command not found.\"" << endl;
      break;
    case fIllegalApplicationState:
      G4cout << "@@ErrResult \"Illegal application state -- command refused\"" << endl;
      break;
    case fParameterOutOfRange:
      G4cout << "@@ErrResult \"Parameter Out of Range.\"" << endl;
      break;
    case fParameterUnreadable:
      G4cout << "@@ErrResult \"Parameter Unreadable.\"" << endl;
      break;
    case fParameterOutOfCandidates:
      G4cout << "@@ErrResult \"Parameter Out of Candidates.\"" << endl;
      break;
    default:
      G4cout << "@@ErrResult \"command refused (" << commandStatus << ")\"" << endl;
    }
  }
}


G4int G4UIGAG::ReceiveG4cout(G4String coutString)
{
  cout << coutString << flush;
  return 0;
}

G4int G4UIGAG::ReceiveG4cerr(G4String cerrString)
{
  cerr << cerrString << flush;
  return 0;
}                                                    

void G4UIGAG::Prompt(G4String aPrompt)
{
  promptCharacter = aPrompt;
}

G4String G4UIGAG::GetCommand()
{
  G4String newCommand;
  G4String nullString;
  while( 1 )
  {
    G4UIcommandTree* tree = UI->GetTree();
    if ( uiMode != terminal_mode ){
      G4cout << "@@PROMPT \"" << promptCharacter << "\"" << endl;
    }
    if ( uiMode != java_mode ){
      G4cout << promptCharacter << "> " << flush;
    }else{
      G4cout << "@@Ready" << endl;
    }
    newCommand.readLine( cin );
    if (!cin.good()) { cin.clear(); newCommand = nullString; iExit=false;break;}
    while( newCommand(newCommand.length()-1) == '_' )
    {
      G4String newLine;
      newCommand.remove(newCommand.length()-1);
      newLine.readLine( cin );
      if (!cin.good()) { cin.clear(); newCommand = nullString; iExit=false;break;}
      newCommand.append(newLine);
    }

    G4String nC = newCommand.strip(G4String::leading);
    // -------------------- nC.toUpper();
    if( nC == "@@GAGmodeJAVA" ) {
      uiMode = java_mode;
      G4cout << endl << "@@Version " << JVersion << endl;
      SendCommandProperties(tree);
      NotifyStateChange();
    }
    else if( nC == "@@GAGmodeTcl" ) {
      uiMode = tcl_mode;
      G4cout << endl << "@@Version " << TVersion << endl;
      SendCommandProperties(tree);
      NotifyStateChange();
    }
    else if( nC(0) == '#' )
      { G4cout << nC << endl; }

    else if( nC == "ls"  || nC(0,3) == "ls " )
    { ListDirectory( nC ); }
    else if( nC == "pwd" )
    { G4cout << "Current Working Directory : " << prefix << endl; }
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
      { G4cout << i << ": " << UI->GetPreviousCommand(i) << endl; }
    }
    else if( nC(0) == '!' )
    {
      G4String ss = nC(1,nC.length()-1);
      G4int vl;
      const char* tt = ss;
      istrstream is((char*)tt);
      is >> vl;
      G4int nh = UI->GetNumberOfHistory();
      if(vl>=0 && vl<nh)
      { 
        newCommand = UI->GetPreviousCommand(vl); 
        G4cout << newCommand << endl;
        break;
      }
      else
      { G4cerr << "history " << vl << " is not found." << endl; }
    }
    else if( nC(0,4) == "exit" )
    { 
      if( iCont )
      { 
	if ( uiMode == terminal_mode){
	  G4cerr << "You are now processing RUN." << endl;
	  G4cerr << "Please abrot it using \"/run/abort\" command first" << endl;
	  G4cerr << " and use \"continue\" command until the application" << endl;
	  G4cerr << " becomes to Idle." << endl;
	}else{
	  G4cout << "@@ErrResult \"You are now processing RUN.\"" << endl;
	}
      }
      else
      {
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

G4String G4UIGAG::GetFullPath( G4String aNewCommand )
{
  G4String newCommand = aNewCommand.strip(G4String::both);
  G4String tmpString;
  if( newCommand(0) == '/' ) 
  { tmpString = newCommand; }
  else if( newCommand(0,3) == "../" )
  {
    G4String tmpPrefix = prefix;
    G4int i_direc = 0;
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

void G4UIGAG::SessionTerminate()
{
  G4cout << "***** Terminal session end *****" << endl;
}

void G4UIGAG::ShowCurrent( G4String newCommand )
{
  G4String theCommand = GetFullPath(newCommand(1,newCommand.length()-1));
  G4String curV = UI->GetCurrentValues(theCommand);
  if( ! (curV.isNull()||curV(0)=='\0' ) ) {
    if (uiMode == terminal_mode){
      G4cout << "Current value(s) of the parameter(s) : " << curV << endl;
    }else{
      G4cout << "@@CurrentValue " << curV << endl;
    }
  } else if (uiMode == terminal_mode){
      G4cout << "Current value is not available." << endl;
    } else {
      G4cout << "@@ErrResult \"Current value is not available.\"" << endl;
    }
}

void G4UIGAG::ChangeDirectory( G4String newCommand )
{
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
    G4cout << "Directory <" << prefix << "> is not found." << endl;
    prefix = savedPrefix;
  }
}

void G4UIGAG::ListDirectory( G4String newCommand )
{
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
  { G4cout << "Directory <" << targetDir << "> is not found." << endl; }
  else
  { commandTree->ListCurrent(); }
}

void G4UIGAG::TerminalHelp(G4String newCommand)
{
  G4UIcommandTree * treeTop = UI->GetTree();
  int i = newCommand.index(" ");
  if( i != RW_NPOS )
  {
    G4String newValue = newCommand(i+1,newCommand.length()-(i+1));
    newValue.strip(G4String::both);
    if( newValue(0) != '/' )
    { newValue.prepend( prefix ); }
    G4UIcommand * theCommand = treeTop->FindPath( newValue );
    if( theCommand != NULL ) 
    { 
      theCommand->List();
      return;
    }
    else
    {
      G4cout << "Command <" << newValue << " is not found." << endl;
      return;
    }
  }

  G4UIcommandTree * floor[10];
  floor[0] = treeTop;
  int iFloor = 0;
  int prefixIndex = 1;
  while( prefixIndex < prefix.length()-1 )
  {
    int ii = prefix.index("/",prefixIndex);
    floor[iFloor+1] =
      floor[iFloor]->GetTree(prefix(0,ii+1));
    prefixIndex = ii+1;
    iFloor++;
  }
  floor[iFloor]->ListCurrentWithNum();
  // 1998 Oct 2 non-number input
  while(1){
    G4cout << endl << "Type the number ( 0:end, -n:n level back ) : "<<flush;
    cin >> i;
    if(!cin.good()){
      cin.clear();
      cin.ignore(30,'\n');
      G4cout << endl << "Not a number, once more" << endl; continue;}
    else if( i < 0 ){
      iFloor += i;
      if( iFloor < 0 ) iFloor = 0;
      floor[iFloor]->ListCurrentWithNum(); continue;}
    else if(i == 0) { break;}
    else if( i > 0 ) {
      int n_tree = floor[iFloor]->GetTreeEntry();
      if( i > n_tree )
      {
        if( i <= n_tree + floor[iFloor]->GetCommandEntry() )
        {
          floor[iFloor]->GetCommand(i-n_tree)->List();
          //iFloor++;
        }
      }
      else
      {
        floor[iFloor+1] = floor[iFloor]->GetTree(i);
        iFloor++;
        floor[iFloor]->ListCurrentWithNum();
      }
    }

  }
  G4cout << "Exit from HELP." << endl << endl;
  G4cout << endl;
  // cin.flush();
  char temp[100];
  cin.getline( temp, 100 );
}







G4String G4UIGAG::ModifyPrefix(G4String newCommand)
{
  G4String newPrefix = prefix;
  while( 1 )
  {
    if( newCommand(0,2) == ".." )
    {
      if( newPrefix != "/" )
      { 
	G4String tmpString = newPrefix(0,newPrefix.length()-1);
        newPrefix = newPrefix(0,tmpString.last('/')+1); 
      }
    }
    else
    {
      newPrefix += newCommand;
      break;
    }
    if( newCommand == ".." || newCommand == "../" )
    { break; }
    newCommand = newCommand(3,newCommand.length()-3);
  }
  return newPrefix;
}

G4UIcommandTree * G4UIGAG::FindDirPath(G4String newCommand)
{
  G4UIcommandTree * comTree = UI->GetTree();
  int idx = 1;
  while( idx < newCommand.length()-1 )
  {
    int i = newCommand.index("/",idx);
    comTree = comTree->GetTree(newCommand(0,i+1));
    if( comTree == NULL ) 
    { return NULL; }
    idx = i+1;
  }
  return comTree;
}

// ----- for JAVA GAG (by T.Kodama)

void G4UIGAG::SendCommandProperties(G4UIcommandTree * tree)
{
  if( tree == NULL ) { 
    G4cerr << "GetTree() returnes null." << endl;
    return;
  }
  if (uiMode == java_mode){
    G4cout << "@@JTreeBegin" << endl;
    CodeGenJavaTree(tree, 0);  
    G4cout << "@@JTreeEnd" << endl;
    CodeGenJavaParams(tree, 0);
  }else{
    G4cout << endl << "@@maketree_start" << endl;
    CodeGenTclTree(tree,0);  
    G4cout << "@@maketree_end" << endl;
    CodeGenTclParams(tree, 0);
  }
}
void G4UIGAG::SendParameterProperties(G4UIcommandTree * tree)
{
  if( tree == NULL ) { 
    G4cerr << "GetTree() returnes null." << endl;
    return;
  }
  if (uiMode == java_mode){
    CodeGenJavaParams(tree, 0);
  }else{
    CodeGenTclParams(tree, 0);
  }
}

void G4UIGAG::CodeGenJavaTree(G4UIcommandTree * tree, int level)
{ 
  int treeEntry, commandEntry;
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  if(level!=0) {
    for(int i=0; i<commandEntry; i++){
      G4cout << tree->GetCommand(i+1)->GetCommandPath() << endl;
    }
  }
  if(treeEntry == 0) return; //end recursion
  
  for(int j=0; j<treeEntry; j++){
    CodeGenJavaTree(tree->GetTree(j+1), level+1);
  }
}

void G4UIGAG::CodeGenJavaParams(G4UIcommandTree * tree, int level) //recursive
{
  int treeEntry, commandEntry, i;
  G4UIcommand * Comp;
  G4UIcommandTree * treeLink;

  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  for(i=0; i<commandEntry; i++) {
    SendAParamProperty(tree->GetCommand(i+1));
  }
  if( treeEntry == 0 )  return;     // end recursion

  for(i=0; i< treeEntry; i++) {
    treeLink = tree->GetTree(i+1);
    G4cout << "@@JDirGuideBegin" << endl;
    G4cout << treeLink->GetPathName() << endl << treeLink->GetTitle() << endl;
    G4cout << "@@JDirGuideEnd" << endl;
    CodeGenJavaParams(treeLink, level+1); 
  }
}

void G4UIGAG::SendAParamProperty(G4UIcommand * Comp)
{
  int guidanceEntry, parameterEntry;
  G4String title, title2;
  G4UIparameter * prp; 
  char c[2];
  guidanceEntry = Comp->GetGuidanceEntries();
  parameterEntry = Comp->GetParameterEntries();
  G4cout << "@@JParamBegin" << endl;
  G4cout << Comp->GetCommandPath() << endl;
  G4cout << guidanceEntry << endl;
  for (int j=0; j<guidanceEntry; j++){
    title = Comp->GetGuidanceLine(j);
    title2 = "";
    if (title != ""){
      for(int i=0; i< title.length(); i++){
	c[0]=title(i);
	c[1]= '\0';
	if ( c[0] == '\n' || c[0] == '\r') {
	  c[0]= ' ';
	}
	title2.append(c);
      }
    }
    G4cout << title2 << endl;
  }
  G4cout << Comp->GetRange() << endl;
  G4cout << parameterEntry << endl;
  for( int par=0; par<parameterEntry; par++) {
    prp = (G4UIparameter *)Comp->GetParameter(par);
    G4cout << prp->GetParameterName() << endl;
    G4cout << prp->GetParameterGuidance() << endl;
    G4cout << prp->GetParameterType() << endl;
    G4cout << prp->IsOmittable() << endl;
    G4cout << prp->GetDefaultValue() << endl;
    G4cout << prp->GetParameterRange() << endl;
    G4cout << prp->GetParameterCandidates() << endl;
  }
  G4cout << "@@JParamEnd" << endl;
}

void G4UIGAG::SendDisableList(G4UIcommandTree * tree, int level)
{ 
  int treeEntry, commandEntry;
  G4UIcommand * Comp;
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  for(int com=0; com<commandEntry; com++) {
    Comp = tree->GetCommand(com+1);
    if( Comp->IsAvailable()==false ) {
       G4cout << Comp->GetCommandPath()<<endl;
    }
  }
  if( treeEntry == 0 ) return;     // end recursion

  for( int i=0; i<treeEntry; i++) {
    SendDisableList(tree->GetTree(i+1), level+1); 
    // be sure the function name is the same
  }
}

// ----- for Tcl GAG

void G4UIGAG::CodeGenTclTree(G4UIcommandTree * tree, int level)
{ 
  int i, j;
  int treeEntry, commandEntry, guidanceEntry;
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();
  G4String commandPath, pathName, title1, title2;
  G4UIcommandTree * t;
  G4UIcommand * Comp;

  for(int com=0; com<commandEntry; com++){
    Comp = tree->GetCommand(com+1);
    commandPath = Comp->GetCommandPath();
    G4cout << commandPath << " @@command" << endl;
    guidanceEntry = Comp->GetGuidanceEntries();
    if (guidanceEntry == 0){
      title2 = "...Title not available...";
    } else {
      title2 = "";
      j = 0;
      while(1){
	title1 = Comp->GetGuidanceLine(j);
	for(i=0; i< title1.length(); i++){
	  char c[2];
	  c[0]=title1(i);
	  c[1]= '\0';
	  if( c[0] == '\"') {
	    title2.append("\\\""); // a Backslash and a double quote
	  } else if ( c[0] == '\n' || c[0] == '\r') {
	    title2.append("\\n");
	  } else title2.append(c);
	}
	j++;
	if (j >= guidanceEntry) break;
	title2.append("\\n");
      }
    }
    G4cout << commandPath << " @@title \""<< title2 <<"\""<< endl;
  }

  if(treeEntry == 0) return; //end recursion

  for(i=0; i< treeEntry; i++){
    t = tree->GetTree(i+1);
    pathName =  t->GetPathName();   
    title1 = t->GetTitle();
    title2 = "";
    for(int i=0; i<title1.length(); i++){
      char c[2];
      c[0]=title1(i);
      c[1]= '\0';
      if( c[0] == '\"') 
	title2.append("\\\""); // a Backslash and a double quote
      else title2.append(c);
    }
    if(level==0) G4cout << pathName<< endl;
    else G4cout << pathName<< "  @@cascade"<<endl;
    G4cout << pathName << "  @@title \"" << title1  << "\""<<endl;
    CodeGenTclTree(t, level+1);
  }
}

void G4UIGAG::CodeGenTclParams( G4UIcommandTree * tree, int level) // recursive
{
  int treeEntry, commandEntry;
  G4UIcommand * Comp;
  treeEntry = tree->GetTreeEntry();
  commandEntry = tree->GetCommandEntry();

  for(int com=0; com<commandEntry; com++) {
    Comp = tree->GetCommand(com+1);
    SendATclParamProperty(Comp);
  }
  if( treeEntry == 0 ) return;     // end recursion

  for( int i=0; i<treeEntry; i++) {
    CodeGenTclParams(tree->GetTree(i+1), level+1); 
    // be sure the function name is the same
  }
}

void G4UIGAG::SendATclParamProperty(G4UIcommand * Comp)
{
    G4UIparameter * prp; 
    int parameterEntry = Comp->GetParameterEntries();
    G4String commandPath = Comp->GetCommandPath();
    G4String commandRange = Comp->GetRange();
    G4cout << "@@parameter_start" << endl;
    G4cout << commandPath << " @@param " << parameterEntry << endl;
    G4cout << "@@command_range \"" << commandRange << "\"" << endl;
    for( int par=0; par<parameterEntry; par++) {
      prp = (G4UIparameter *)Comp->GetParameter(par);
      G4cout << "{" ;
      G4cout << "@@param_name : \"" << prp->GetParameterName() <<"\""<<endl;
      G4String  guide1,guide2;
      guide1 = prp->GetParameterGuidance();
      guide2 = "";
      for(int i=0; i<guide1.length(); i++){
        char c[2];
        c[0]=guide1(i);
        c[1]= '\0';
        if( c[0] == '\"') 
        guide2.append("\\\""); // a Backslash and a double quote
        else guide2.append(c);
      }
      G4cout << " @@param_guide : \"" << guide2 << "\""<<endl; 
      G4cout << " @@param_type : \"" << prp->GetParameterType()<<"\""<<endl;
      G4cout << " @@param_omit : \"" << prp->IsOmittable()<<"\""<<endl;
      G4cout << " @@param_default : \""<< prp->GetDefaultValue()<<"\""<<endl;
      G4cout << " @@param_range : \""<< prp->GetParameterRange()<<"\""<<endl;
      G4cout << " @@param_candidate : \"" << prp->GetParameterCandidates()<< "\""<<endl;
      G4cout << "}" << endl;
    }
    G4cout << "@@parameter_end" << endl;
}

void G4UIGAG::NotifyStateChange(void)
{
   G4String stateString;
   G4StateManager * statM = G4StateManager::GetStateManager();
   G4UIcommandTree * tree = UI->GetTree();
   stateString = statM->GetStateString(statM->GetCurrentState());
   if ( uiMode != terminal_mode ){
     G4cout << "@@State \"" << stateString << "\"" << endl;
     G4cout << "@@DisableListBegin"<<endl;
     SendDisableList(tree, 0);
     G4cout << "@@DisableListEnd" <<endl;
   }
}

void G4UIGAG::NotifyCommandUpdate(void)
{
  G4UIcommandTree * tree = UI->GetTree();
  SendCommandProperties(tree);
}

void G4UIGAG::NotifyParameterUpdate(G4UIcommand* com)
{
  if (uiMode == java_mode) 
    SendAParamProperty(com);
  else
    SendATclParamProperty(com);
}

//####### update check routines ####################################
void G4UIGAG::UpdateState(void)
{
   static G4ApplicationState previousState= PreInit;
   G4ApplicationState  newState;
   G4StateManager *statM = G4StateManager::GetStateManager();
   newState = statM->GetCurrentState();
   if( newState != previousState ) 
   {
      NotifyStateChange();
      previousState = newState; 
   }
}

int G4UIGAG::CommandUpdated(void)
{
  int added=0, deleted=0;
  int pEntry= previousTreeCommands.entries();
  int nEntry= newTreeCommands.entries();
  int i,j;
  for( i=0; i<pEntry; i++) {      // check deleted command(s)
      for( j=0; j<nEntry; j++) {
         if( previousTreeCommands(i) == newTreeCommands(j)) break;
      }
      if( j==nEntry ) { 
         deleted = 1;
         //G4cout <<"deleted: "<< previousTreeCommands(i) << endl;
      }
  }
  for( i=0; i<nEntry; i++) {      // check added command(s)
      for( j=0; j<pEntry; j++) {
         if( newTreeCommands(i) == previousTreeCommands(j)) break;
      }
      if( j==pEntry ) { 
         added = 1;
      //   G4cout <<"added: "<< newTreeCommands(i) << endl;
      }
  }
  if( added    && deleted==0 ) {G4cout<<"c added"<<endl;return added;}
  if( added==0 && deleted ) {G4cout<<"c deleted"<<endl;return deleted;}
  if( added    && deleted ) {G4cout<<"c add/deleted"<<endl;return addedAndDeleted;}
  return notChanged;
}

void G4UIGAG::GetNewTreeStructure(G4UIcommandTree * tree, int level)
{ 
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
      newTreeCommands.insert( commandPath + " " + title );
  }

  if(treeEntry == 0) return; //end recursion

  for(int i=0; i< treeEntry; i++){
    t = tree->GetTree(i+1);
    pathName =  t->GetPathName();   
    title = t->GetTitle();
    newTreeCommands.insert( pathName + " " + title );
    GetNewTreeStructure(t, level+1);
  }
}

void G4UIGAG::UpdateParamVal(void)
{
  // call NotifyParameterUpdate() if the value of each
  //  command/parameter is updated.
  //  assuming the command structure is not changed.
  int pEntry= previousTreeParams.entries();
  int nEntry= newTreeParams.entries();
  int i;
  G4UIcommand* Comp;
  if (pEntry != nEntry) return; 
  for( i=0; i<nEntry; i++) {
    if( previousTreeParams(i) != newTreeParams(i)){
       Comp = newTreePCP(i);
       G4cout << Comp->GetCommandPath()
            << " command is updated." <<endl; 
       NotifyParameterUpdate(Comp);
    }
  }
}

//void G4UIGAG::paramUpdate(void)
//{  
//  int added=0, deleted=0;
//  int pEntry= previousTreeParams.entries();
//  int nEntry= newTreeParams.entries();
//  int i,j;
//
//  if (pEntry != nEntry)  return NULL;
//  for( i=0; i<pEntry; i++) {       // check deleted param(s)
//    for( j=0; j<nEntry; j++) {
//       if( previousTreeParams(i) == newTreeParams(j)) break;
//    }
//    if( j==nEntry ) { 
//       deleted = 1;
//       //G4cout <<"para deleted: "<< previousTreeParams(i) << endl;
//    }
//  }
//  for( i=0; i<nEntry; i++) {      // check added param(s)
//    for( j=0; j<pEntry; j++) {
//       if( newTreeParams(i) == previousTreeParams(j)) break;
//    }
//    if( j==pEntry ) { 
//       added = 1;
//       //G4cout <<"para added: "<< newTreeParams(i) << endl;
//    }
//  }
//  if( added    && deleted==0 ) {G4cout<<"p added"<<endl;return added;}
// if( added==0 && deleted )  {G4cout<<"p deleted"<<endl;return deleted;}
//  if( added    && deleted )  {G4cout<<"p add/deleted"<<endl; return addedAndDeleted;}
//  return notChanged;
//}

void G4UIGAG::GetNewTreeValues( G4UIcommandTree * tree, int level) // recursive
{
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
     newTreeParams.insert( param + "\n"); 
     newTreePCP.insert( Comp ); 
   }
   if( treeEntry == 0 )  return;     // end recursion
   for( int i=0; i< treeEntry; i++) {
      t = tree->GetTree(i+1);
      GetNewTreeValues(t, level+1);
   }
}
//######################################################
