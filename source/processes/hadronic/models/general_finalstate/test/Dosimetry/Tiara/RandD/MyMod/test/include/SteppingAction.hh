#ifndef __STEPPING_ACTION_DEFINED__
#define __STEPPING_ACTION_DEFINED__

#include "G4UserSteppingAction.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "stdio.h"

class G4Step;
class Analyzer;
class Hall;
class RunAction;
class G4Event;
class EventAction;


class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(Analyzer* pAnalyzer,Hall* pDetector,RunAction* pRunAction);
  ~SteppingAction();
  void UserSteppingAction(const G4Step* pStep);
  void SetEvtAction(EventAction* pEAction){m_pEvtAction=pEAction;};
  void ForcedBeginEvent();
  void ForcedEndEvent();
  void UseFile(char* szFileName);
  G4double GetSigmaFromFile();
private:
  Analyzer* m_pAnalyzer;
  Hall*     m_pDetector;
  RunAction* m_pRunAction;
  EventAction* m_pEvtAction;
  double m_dNi;
  FILE* m_pFile;
  char* m_szFileName;
};

class EventAction : public G4UserEventAction
{
  unsigned m_nCurrID;
  SteppingAction* m_pStepAction;
public:
  EventAction():m_nCurrID(0),m_pStepAction(NULL){};
  void BeginOfEventAction(const G4Event* pEvt);
  void EndOfEventAction(const G4Event* pEvt);
  void ResetCounter(){m_nCurrID=0;};
  unsigned GetCounter(){return m_nCurrID;};
  void SetStepAction(SteppingAction* pStepAction){m_pStepAction = pStepAction;};
};
#endif
