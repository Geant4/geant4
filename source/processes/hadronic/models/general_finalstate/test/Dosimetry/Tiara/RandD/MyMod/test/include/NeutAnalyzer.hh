#ifndef __ANALYZER_DEFINED__
#define __ANALYZER_DEFINED__

#define MAX_HISTOS 13

#include "globals.hh"

enum {enHistoAtBeam=1,enHisto20Cm=2,enHisto40Cm=4,enHistoGhost=8};
enum {enEnergy43=1,enEnergy68=2};

class AnalyzerMessenger;
class RunAction;
class ParticleGun;

class Analyzer
{
public:
  Analyzer(RunAction* pAction,ParticleGun* pGun,G4bool bShow=true);
  ~Analyzer();
  void Histograms(unsigned char eNumerator);
  void Fill(unsigned char What,G4double value,G4double weight=1.,G4double EvWeight=1.);
  void BuildDifferences(unsigned char What);
  void Reset();
  int GetAfterTargetID(){ return HistoArray[3];};
  void ShowHisto(unsigned char What);
  void HideHisto(unsigned char What);
  void MakePS(unsigned char What,G4String szFileName);
  void SaveHisto(unsigned char What);
  void DeleteHisto(unsigned char What);
  void AddData(unsigned char What,G4String szValue);
  void AdditionalHisto(unsigned char What,G4String storeName,unsigned StoreID);
  void SetProperty(unsigned char What,G4String szKey,G4String szValue);
  void ReadOptions(unsigned char What,G4String szFileName);
  void DeleteOption(unsigned char What,G4String szKey);
  void ListOptions(unsigned char What);
  void UseCut(double dCut){m_dCutEnergy = dCut;};
  void SetMinEnergy(double dMinEnergy){m_dMinEnergy = dMinEnergy;};
  void BuildSmooth(unsigned char What);
  void Rebin(unsigned char First,unsigned char Last,unsigned char bias);
  double GetEnergy(unsigned char Which,double Energy);
  double GetEnergyRange(unsigned char Which,double Energy);
private:
  double        m_dCutEnergy;
  double        m_dMinEnergy;
  unsigned char m_En;
  G4int         HistoArray[MAX_HISTOS];
  RunAction*    m_pAction;
  AnalyzerMessenger* m_pMessenger;
  ParticleGun*  m_pGun;
};
#endif
