#include "RunAction.hh"
#include "G4Run.hh"
#include "AnalyzerLin.hh"
#include "NeutAnalyzer.hh"
#include "Hall.hh"

//#include "stdafx.h"

RunAction::RunAction(ParticleGun* pGun,Hall* pHall)
{
  m_pHall = pHall;
  m_pGun = pGun;
  m_cProtons = 0;
  m_cOldProtons = 1.;
}

RunAction::~RunAction()
{
}
void RunAction::BeginOfRunAction(const G4Run* pRun)
{
}
void RunAction::EndOfRunAction(const G4Run* pRun)
{
  if(m_pGun->GetStat())
    m_cProtons = alCalculateIntegral(m_pAnalyzer->GetAfterTargetID())*
      m_cOldProtons/(m_pGun->GetRArea());
  else
    m_cProtons += eplus*1e-6*coulomb*pRun->GetNumberOfEvent();
  alScaleAllHistos((G4double)m_cOldProtons/(G4double)m_cProtons);
  m_cOldProtons = m_cProtons;
}
void RunAction::ResetCounter()
{
  m_cProtons = 0;
  m_cOldProtons = 1.;
}

G4double RunAction::GetNeutrons()
{
  G4double temp;
  temp = alGetMinimumEntries(m_pAnalyzer->GetAfterTargetID()-1);
  G4cout<<"Minimum entries: "<<temp<<G4endl;
  return temp;
}
G4double RunAction::GetAxisNeutrons()
{
  return alGetMinimumEntries(m_pAnalyzer->GetAfterTargetID()-3);
}
