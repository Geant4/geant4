#ifndef __EVENT_ACTION_DEFINED__
#define __EVENT_ACTION_DEFINED__

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "ParticleGun.hh"

class Analyzer;
class G4Run;
class Hall;

class RunAction : public G4UserRunAction
{
public:
  RunAction(ParticleGun* pGun,Hall* pHall);
  ~RunAction();
  void BeginOfRunAction(const G4Run* pEvent);
  void EndOfRunAction(const G4Run* pEvent);
  void ResetCounter();
  G4double GetProtonsNumber(){
    return (m_cOldProtons);
  };
  void Register(Analyzer* pAnalyzer){m_pAnalyzer = pAnalyzer;};
  G4double GetNeutrons();
  G4double GetAxisNeutrons();
  bool GetStat(){return m_pGun->GetStat();};
private:
  ParticleGun* m_pGun;
  Hall*        m_pHall;
  G4double m_cProtons;
  G4double m_cOldProtons;
  Analyzer* m_pAnalyzer;
};

#endif
