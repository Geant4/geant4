#include "NeutAnalyzer.hh"
#include "AnalyzerLin.hh"
#include "AnalyzerMessenger.hh"
#include "RunAction.hh"
#include "ParticleGun.hh"

//#include "stdafx.h"

Analyzer::Analyzer(RunAction* pAction,ParticleGun* pGun,G4bool bShow)
{
  m_pGun = pGun;
  if(m_pGun->GetStat())
    m_En = ((m_pGun->GetMax() < 50.) ? enEnergy43 : enEnergy68);
  else
    m_En = (m_pGun->CurrEnergy() < 50.) ? enEnergy43 : enEnergy68;
  m_dCutEnergy = 40;
  m_dMinEnergy = 10;
  m_pAction = pAction;
  m_pAction->Register(this);
  for(int i=0;i<MAX_HISTOS;i++)
    HistoArray[i]=-1;
  alInit();
  m_pMessenger = new AnalyzerMessenger(this);
  ShowViewer(bShow);
}

Analyzer::~Analyzer()
{
  alClear();
  delete m_pMessenger;
}

void Analyzer::Histograms(unsigned char eNumerator)
{
  int index=0;
  char title[512],*help;
  int nBins=100;
  char* xTitle,*yTitle;
  double min,max;
  if(m_pGun->GetStat())
    m_En = (m_pGun->GetMax() < 50.) ? enEnergy43 : enEnergy68;
  else
    m_En = (m_pGun->CurrEnergy() < 50.) ? enEnergy43 : enEnergy68;
  while(index<4){
    if(HistoArray[index]!=-1){
      alDeleteHisto(HistoArray[index]);
      HistoArray[index]=-1;
    }
    if(eNumerator& ~(1<<index)){
      switch(index){
      case 0: help = " at beam axis"; min = 3e+6; max = 1e+8; 
	xTitle = "Energy (eV)"; 
	yTitle = "Flux (n sr^-1! [m]C^-1!leth^-1!)";
	break;
      case 1: help = " 20 cm from beam axis"; min = 3e+6; max = 1e+8;
	xTitle = "Energy (eV)";
	yTitle = "Flux (n sr^-1! [m]C^-1!leth^-1!)";
	break;
      case 2: help = " 40 cm from beam axis"; min= 3e+6; max = 1e+8;
	xTitle = "Energy (eV)";
	yTitle = "Flux (n sr^-1! [m]C^-1!leth^-1!)";
	break;
      case 3: help = " after target"; min = 0; max = 80;
	xTitle = "Energy (MeV)";
	yTitle = "Flux (n sr^-1! [m]C^-1!MeV^-1!)";
	break;
      }
      sprintf(title,"For protons at %2d MeV",(m_En==enEnergy43) ? 43 : 68);
      alSetTitle(title);
      sprintf(title,"Neutron Flux %s",help);
      HistoArray[index] = alAddHisto(title,nBins,min,max);
      alAddTitles(HistoArray[index],xTitle,yTitle);
    }
    index++;
  }
}
void Analyzer::Fill(unsigned char What,G4double value,G4double weight,G4double EvWeight)
{
  switch(What){
  case 1: What=0; break;
  case 2: What = 1; break;
  case 4: What = 2; break;
  case 8: What = 3; break;
  }
  if(HistoArray[What]!=-1){
    alFill(HistoArray[What],value,weight,EvWeight);
    alAddEvent(HistoArray[What]);
  }
}
void Analyzer::BuildDifferences(unsigned char What)
{
  char title[512],*help;
  G4int nBins = 100;
  G4double min,max,dCut,dCut1;
  if(m_pGun->GetStat())
    m_En = (m_pGun->GetMax() < 50) ? enEnergy43 : enEnergy68;
  else
    m_En = (m_pGun->CurrEnergy() < 50) ? enEnergy43 : enEnergy68;
  switch(What){
  case 1: What = 0;
    help = "beam axis";
    min = 3e+6;
    max = 1e+8;
    dCut = m_dCutEnergy/eV;
    dCut1 = m_dMinEnergy/eV;
    break;
  case 2: What = 1;
    help = "20 cm from beam axis";
    min = 3e+6; max = 1e+8;
    dCut1 = m_dMinEnergy/eV;
    dCut = m_dCutEnergy/eV;
    break;
  case 4: What = 2;
    help = "40 cm from beam axis";
    min = 3e+6;
    max = 1e+8;
    dCut = m_dCutEnergy/eV;
    dCut1 = m_dMinEnergy/eV;
    break;
  case 8: What = 3;
    help = "after the target";
    min = 0.;
    max = 80.;
    dCut = m_dCutEnergy;
    dCut1 = m_dMinEnergy;
    break;
  }
  if(HistoArray[What]==-1) return;
  if(HistoArray[What+4]==-1){
    sprintf(title,"Differences(relative) %s",help);
    HistoArray[What+4] = alAddHisto(title,nBins,min,max);
  }
  alSubstract(HistoArray[What],HistoArray[What+4],dCut,dCut1);
}

void Analyzer::Reset()
{
  m_pAction->ResetCounter();
  alReset();
  alRefreshHisto(false);
}

void Analyzer::ShowHisto(unsigned char What)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alShowHisto(HistoArray[What-1],0);
}

void Analyzer::HideHisto(unsigned char What)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alHideHisto(HistoArray[What-1]);
}

void Analyzer::MakePS(unsigned char What,G4String szFileName)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alMakePS(HistoArray[What-1],szFileName);
}

void Analyzer::SaveHisto(unsigned char What)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alSaveHisto(HistoArray[What-1]);
}

void Analyzer::AddData(unsigned char What,G4String szValue)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alAddDataFromFile(HistoArray[What-1],szValue);
  alRefreshHisto(false);
}

void Analyzer::DeleteHisto(unsigned char What)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alDeleteHisto(HistoArray[What-1]);
  HistoArray[What-1] = -1;
}

void Analyzer::SetProperty(unsigned char What,G4String szKey,G4String szValue)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID:"<<(unsigned)What<<G4endl;
    return;
  }
  alSetProperty(HistoArray[What-1],szKey,szValue);
  alRefreshHisto(false);
}

void Analyzer::ReadOptions(unsigned char What,G4String FileName)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alReadOptions(HistoArray[What-1],(char*)FileName.data());
}

void Analyzer::DeleteOption(unsigned char What,G4String szKey)
{
  if((What > MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alDeleteOptions(HistoArray[What-1],szKey);
}

void Analyzer::ListOptions(unsigned char What)
{
  if((What> MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  alListOptions(HistoArray[What-1]);
}

void Analyzer::BuildSmooth(unsigned char What)
{
  if((What>MAX_HISTOS)||(HistoArray[What-1]==-1)){
    G4cout<<"Analyzer:: Wrong histo ID: "<<(unsigned)What<<G4endl;
    return;
  }
  if(What>4){
    G4cout<<"Only \"data\" histograms can be smoothed :-)"<<G4endl;
    return;
  }
  unsigned char targetHisto = What+7;
  What--;
  double dCut;
  switch(What){
  case 0:
  case 1:
  case 2: dCut = m_dCutEnergy*1e+6; break;
  default : dCut = m_dCutEnergy;
  }
  if(HistoArray[targetHisto]==-1){
    G4String HistoName = "Smoothed histo ";
    G4String xTitle,yTitle;
    double nBins = 100;
    double dMin,dMax,dCut;
    switch(What){
    case 0: HistoName = HistoName+" at beam axis"; dMin = 3e+6; dMax = 1e+8; 
      xTitle = "Energy (eV)"; 
      yTitle = "Flux (n sr^-1! [m]C^-1!leth^-1!)";
      dCut = m_dCutEnergy*1e+6;
      break;
    case 1: HistoName = HistoName+" 20 cm from beam axis"; 
      dMin = 3e+6; dMax = 1e+8;
      xTitle = "Energy (eV)";
      yTitle = "Flux (n sr^-1! [m]C^-1!leth^-1!)";
      dCut = m_dCutEnergy*1e+6;
      break;
    case 2: HistoName = HistoName+" 40 cm from beam axis"; 
      dMin= 3e+6; dMax = 1e+8;
      xTitle = "Energy (eV)";
      yTitle = "Flux (n sr^-1! [m]C^-1!leth^-1!)";
      dCut = m_dCutEnergy*1e+6;
      break;
    case 3: HistoName = HistoName+" after target"; 
      dMin = 0; dMax = 80;
      xTitle = "Energy (MeV)";
      yTitle = "Flux (n sr^-1! [m]C^-1!MeV^-1!)";
      dCut = m_dCutEnergy;
      break;
    }
    HistoArray[targetHisto]=alAddHisto(HistoName,(int)(nBins+0.5),dMin,dMax);
    alAddTitles(HistoArray[targetHisto],(char*)xTitle.data(),(char*)yTitle.data());
  }
  alCreateSmooth(HistoArray[What],HistoArray[targetHisto],dCut);
}

void Analyzer::AdditionalHisto(unsigned char What,G4String storeName,unsigned StoreID)
{
  char* name = (storeName=="") ? NULL:storeName.data();
  char label[32];
  label[0]=0;
  int lastChar=0;
  while(StoreID!=0){
    char c = StoreID%10;
    c += '0';
    label[lastChar] = c;
    label[lastChar+1] = 0;
    lastChar++;
    StoreID /= 10;
  }
  lastChar--;
  for(int i=0;i<lastChar;i++,lastChar--){
    label[lastChar] ^= label[i];
    label[i] ^= label[lastChar];
    label[lastChar] ^= label[i];
  }
  if(HistoArray[What-1]!=-1){
    alAdditionalHisto(HistoArray[What-1],name,label);
  }
}

void Analyzer::Rebin(unsigned char First,unsigned char Last,unsigned char bias)
{
  First--;
  Last--;
  unsigned i;
  for(i=0;i<(Last-First+1);i++){
    if((HistoArray[First+i]!=-1) && (HistoArray[First+bias+i]!=-1))
      alRebin(HistoArray[First+i],HistoArray[First+i],HistoArray[First+bias+i]-HistoArray[First+i]);
  }
}

double Analyzer::GetEnergy(unsigned char Which,double Energy)
{
  Which--;
  if(HistoArray[Which]==-1) return 0;
  return alGetEnergy(HistoArray[Which],Energy);
}
double Analyzer::GetEnergyRange(unsigned char Which,double Energy)
{
  Which--;
  if(HistoArray[Which]==-1) return 0;
  return alGetEnergyRange(HistoArray[Which],Energy);
}
