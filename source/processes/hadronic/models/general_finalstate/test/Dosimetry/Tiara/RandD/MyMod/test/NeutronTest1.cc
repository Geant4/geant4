#include "stdafx.h"

RunAction* g_pRunAction;
bool g_pInterrupt = false;
Analyzer* g_pAnalyzer;

class RunManager : public G4RunManager
{
public:
  RunManager() : G4RunManager() {;};
  ~RunManager() {;};
  void ProcessTill(G4int nProtons,G4int nNeutrons);
};
#include <sys/times.h>
#include <time.h>
#include <stdio.h>

FILE* g_pFile;

void RunManager::ProcessTill(G4int nProtons,G4int nNeutrons)
{
  G4RunManager* pManager = G4RunManager::GetRunManager();
  struct tms time;
  clock_t cl_time = 0;
  clock_t tmp_time;
  unsigned ncycles=0;
  G4double neutrons=0,Neutrons = 0;
  static G4double OldNeutrons = 0;
  Neutrons = OldNeutrons;
  while(/*neutrons < (G4double)nNeutrons*/((double)cl_time)/((double)CLOCKS_PER_SEC)<0.01){
    ncycles++;
    times(&time);
    tmp_time = time.tms_utime;
    pManager->BeamOn(nProtons);
    times(&time);
    cl_time += time.tms_utime - tmp_time;
    G4cout<<"next time: "<<((double)cl_time)/((double)CLOCKS_PER_SEC);
    neutrons = g_pRunAction->GetNeutrons();
    Neutrons += nProtons;
    if(g_pInterrupt){
      g_pInterrupt = false;
      break;
    }
    g_pAnalyzer->BuildDifferences(1);
    g_pAnalyzer->BuildDifferences(2);
    g_pAnalyzer->BuildDifferences(4);
    g_pAnalyzer->BuildDifferences(8);
    alRefreshHisto(false);
    G4double sigma = ((SteppingAction*)GetUserSteppingAction())->GetSigmaFromFile();
    fprintf(g_pFile,"%g\t%g\t%g\n",(time.tms_utime-tmp_time)/(double)CLOCKS_PER_SEC,Neutrons-OldNeutrons,(Neutrons-OldNeutrons)*(double)CLOCKS_PER_SEC/(double)cl_time/sigma);
  }
  fprintf(g_pFile,"\n");
  OldNeutrons = Neutrons;
  fprintf(g_pFile,"ncycles: %d\n",ncycles);
}

RunManager* g_pRunManager;

void TerminateRun(int nSig)
{
  g_pInterrupt = true;
}


main(int argc, char** argv)
{
  FILE* pTmpFile;
  signal(SIGTERM,TerminateRun);
  signal(SIGINT,TerminateRun);
  signal(SIGHUP,SIG_IGN);
  srand(time(NULL)%100+1);
  alInit();
  ParticleGun* pGun = new ParticleGun;
  Hall* pHall = new Hall;
  g_pRunAction = new RunAction(pGun,pHall);
  g_pAnalyzer = new Analyzer(g_pRunAction,pGun);
  g_pRunManager = new RunManager;
  g_pRunManager->SetUserAction(pGun);
  g_pRunManager->SetUserInitialization(new PhysicsList(pHall));
  g_pRunManager->SetUserInitialization(pHall);

  g_pRunManager->SetUserAction(g_pRunAction);
  SteppingAction* pStepAc;
  g_pRunManager->SetUserAction(pStepAc = new SteppingAction(g_pAnalyzer,pHall,g_pRunAction));
  EventAction* pEvAc;
  g_pRunManager->SetUserAction(pEvAc = new EventAction);
  pStepAc->SetEvtAction(pEvAc);
  pEvAc->SetStepAction(pStepAc);
  pStepAc->UseFile(tmpnam(NULL));
  g_pFile = fopen("tmp.txt","wt");

  G4VisManager* pVisManager = new VisManager;
  pVisManager->Initialize();
  
  if(argc==1){
    pTmpFile = fopen("preinit.mac","r");
    if(pTmpFile!=NULL){
      G4String initCmd = G4String("/control/execute preinit.mac");
      G4UImanager::GetUIpointer()->ApplyCommand(initCmd);
      fclose(pTmpFile);
    }
  }
  g_pRunManager->Initialize();

  if(argc > 1){
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    for(G4int i=1;i<argc;i++){
      G4String cmd = G4String("/control/execute ") + G4String(argv[1]);
      UImanager->ApplyCommand(cmd);
    }
  }
  else{
    pTmpFile = fopen("post_init.mac","r");
    if(pTmpFile!=NULL){
      G4String postInitCmd = G4String("/control/execute post_init.mac");
      G4UImanager::GetUIpointer()->ApplyCommand(postInitCmd);
      fclose(pTmpFile);
    }
    G4UIsession* pSession = new G4UIterminal;
    pSession->SessionStart();
    delete pSession;
  }
  //tuk. Pyrwo se izstrivat scorer-ite
  /*  delete pHall;
      g_pRunManager->SetUserInitialization((G4VUserPhysicsList*)NULL);*/
  pHall->DontUseImportance();
  delete pVisManager;
  delete g_pRunManager;
  delete g_pAnalyzer;

  signal(SIGTERM,SIG_DFL);
  signal(SIGINT,SIG_DFL);
}
