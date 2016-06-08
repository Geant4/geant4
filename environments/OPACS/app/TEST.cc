// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TEST.cc,v 1.3 1999/12/15 14:48:40 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
#include <stdlib.h>

/*G4*/
#include <G4RunManager.hh>

/*OPACS*/
/*begin Want_h*/
#include <WoWo.h>
#include <WoXtw.h>
#include <WoXm.h>
#include <WoXo.h>
#include <WoHTML.h>
#include <WoXz.h>
#include <WoDeX.h>
#include <WoCi.h>
#include <WoPacksCi.h>
#include <Wotcl.h>
#include <WoHoXz.h>
#include <WoXaw.h>
#include <WoG3o.h>
#include <WoG4o.h>
#include <WoFNAL.h>
#include <WoPAW.h> /*Need HAS_XZ and HAS_PAW flags.*/
/*end Want_h*/

#include <OHistogram.h>

/*G4o*/
#include <G4oState.hh>

/***************************************************************************/
/***** G4 customization ****************************************************/
/***************************************************************************/
#include <test19DetectorConstruction.hh>
#include <MyPhysicsList.hh>
#include <MyPrimaryGeneratorAction.hh>

#include <G4UserRunAction.hh>
#include <G4Run.hh>
#include <G4Timer.hh>
class RunAction : public G4UserRunAction {
private:
  G4Timer* timer;
  G4int    runIDcounter;
#ifdef HAS_HO
  static   OHistogram myHisto1D;
  static   OHistogram myHisto2D;
#endif /*HAS_HO*/
public:    
  RunAction() {
    timer        = new G4Timer;
    runIDcounter = 0;
  }
  ~RunAction() {
    delete timer;
  }
  void BeginOfRunAction(G4Run* aRun) {
    aRun->SetRunID(runIDcounter++);
    CWarnF ("### Run %d start.\n",aRun->GetRunID());
    timer->Start();
#ifdef HAS_HO
    // Histo are created in an osh script somewhere (G4oRunBegin.osh).
    myHisto1D = OHistogramGetIdentifier("HepJamesRandom1");
    myHisto2D = OHistogramGetIdentifier("HepJamesVsDRand48");
#endif /*HAS_HO*/
  }
  void EndOfRunAction(G4Run* aRun) {
    timer->Stop();
    if (timer->IsValid()) {
      CWarnF ("number of event = %d User=%gs Real=%gs Sys=%gs\n",
	      aRun->GetNumberOfEvent(),
	      timer->GetUserElapsed(),
	      timer->GetRealElapsed(),
	      timer->GetSystemElapsed());
    } else {
      CWarnF ("number of event = %d User=****s Real=****s Sys=****s\n",
	      aRun->GetNumberOfEvent());
    }
#ifdef HAS_HO
    myHisto1D = NULL;
    myHisto2D = NULL;
#endif /*HAS_HO*/
  }
#ifdef HAS_HO
  static OHistogram get_1d(){return myHisto1D; }
  static OHistogram Get2d(){return myHisto2D; }
#endif /*HAS_HO*/
};
#ifdef HAS_HO
OHistogram RunAction::myHisto1D = NULL;
OHistogram RunAction::myHisto2D = NULL;
#endif /*HAS_HO*/

#include <G4UserEventAction.hh>
class EventAction : public G4UserEventAction {
public:
  EventAction() {
  }
  ~EventAction() {
  }
  void BeginOfEventAction() {
  }
  void EndOfEventAction() {
  }
};

#include "G4UserSteppingAction.hh"
class SteppingAction : public G4UserSteppingAction {
private:
  HepJamesRandom theJamesEngine;
  DRand48Engine theDRand48Engine;
public:
  SteppingAction(){};
  ~SteppingAction(){};
  void UserSteppingAction() {
#ifdef HAS_HO
    OHistogram h1 = RunAction::get_1d();
    OHistogram h2 = RunAction::Get2d();
    
    HepRandom::setTheEngine (&theJamesEngine);
    double     james = RandGauss::shoot(0.3,0.1);
    OHistogramFillOneDimensional(h1,james,0.01);

    HepRandom::setTheEngine(&theDRand48Engine);
    double     d48 = RandGauss::shoot(0.7,0.1);
    OHistogramFillTwoDimensional(h2,james,d48,0.01);   
#endif /*HAS_HO*/
  }
};

static char what[] = "@(#)TEST";
/***************************************************************************/
int main (
 int a_argn
,char** a_args
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  what[0] = '@'; /*c++ no warning.*/
/*begin Want_c*/
#include <WoWo.ic>
#include <WoXtw.ic>
#include <WoXm.ic>
#include <WoXo.ic>
#include <WoHTML.ic>
#include <WoXz.ic>
#include <WoDeX.ic>
#include <WoCi.ic>
#include <WoPacksCi.ic>
#include <Wotcl.ic>
#include <WoHoXz.ic>
#include <WoXaw.ic>
#include <WoG3o.ic>
#include <WoG4o.ic>
#include <WoFNAL.ic>
#include <WoPAW.ic>
/*end Want_c*/

  G4RunManager*   runManager = new G4RunManager;
  runManager->SetUserInitialization (new test19DetectorConstruction);
  runManager->SetUserInitialization (new MyPhysicsList); 
  runManager->SetUserAction         (new MyPrimaryGeneratorAction);
  runManager->SetUserAction         (new RunAction);
  runManager->SetUserAction         (new EventAction);
  runManager->SetUserAction         (new SteppingAction);

  G4oState*       state = new G4oState("TEST_");

  WoStartup       (a_argn,a_args);
  WoProcessEvents ();
  WoClearClass    ();

  delete          runManager;
  delete          state;

  return          EXIT_SUCCESS;
}


