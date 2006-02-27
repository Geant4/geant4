// $Id: MyApplication.cc,v 1.1 2006-02-27 10:05:24 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MyApplication.cc
//
//                                         2005 Q
// ====================================================================
#include "MyApplication.hh"
#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "QMaterials.hh"
#include "QDetectorConstruction.hh"
#include "QPhysicsList.hh"
#include "QEventAction.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////
MyApplication::MyApplication()
//////////////////////////////
{
  if(G4RunManager::GetRunManager() ==0 )
    G4RunManager* runManager= new G4RunManager();
}


////////////////////////////////////////////////////////////////////////
MyApplication::MyApplication(int argc, char** argv, const G4String& ver)
////////////////////////////////////////////////////////////////////////
{
  ParseArguments(argc, argv);

  if(G4RunManager::GetRunManager() ==0 )
    G4RunManager* runManager= new G4RunManager();
}


///////////////////////////////
MyApplication::~MyApplication()
///////////////////////////////
{
}


/////////////////////////////////////////////////////////
void MyApplication::ParseArguments(int argc, char** argv)
/////////////////////////////////////////////////////////
{
  if (argc!=1) {
    batchname= argv[1];
  }
}


////////////////////////////////////
void MyApplication::ShowHelp() const
////////////////////////////////////
{
  G4cout << version << G4endl;
}


//////////////////////////////////
void MyApplication::StartSession()
//////////////////////////////////
{
  G4UItcsh* tcsh= new
    G4UItcsh("[40;01;33mq00[40;31m(%s)[40;36m[%/][00;01;30m:");
  tcsh-> SetLsColor(BLUE,RED);

  G4UIterminal* session= new G4UIterminal(tcsh);
  session-> SessionStart();
  delete session;
}


//////////////////////////////////
void MyApplication::ExecuteBatch()
//////////////////////////////////
{
  G4UImanager* UImanager= G4UImanager::GetUIpointer();
  G4String command = "/control/execute ";
  UImanager-> ApplyCommand(command+batchname);
}


///////////////////////////////
void MyApplication::Configure()
///////////////////////////////
{
  G4RunManager* runManager= G4RunManager::GetRunManager();

  // define materials
  QMaterials* materialStore= new QMaterials();
  materialStore-> Construct();

  runManager-> SetUserInitialization(new QDetectorConstruction);
  runManager-> SetUserInitialization(new QPhysicsList);
  //runManager-> SetUserAction(new QEventAction);

  runManager-> Initialize();
}

