#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "NeutAnalyzer.hh"
#include "Hall.hh"
#include "RunAction.hh"

//#include "stdafx.h"
#include "G4Event.hh"

void EventAction::BeginOfEventAction(const G4Event* pEvt)
{
  m_nCurrID++;
  if(m_pStepAction!=NULL) m_pStepAction->ForcedBeginEvent();
}

void EventAction::EndOfEventAction(const G4Event* pEvt)
{
  if(m_pStepAction!=NULL) m_pStepAction->ForcedEndEvent();
}

SteppingAction::SteppingAction(Analyzer* pAnalyzer,Hall* pDetector,RunAction* pRunAction)
{
  m_pAnalyzer = pAnalyzer;
  m_pDetector = pDetector;
  m_pRunAction = pRunAction;
  m_pEvtAction = NULL;
  m_dNi = 0;
  m_pFile = NULL;
  m_szFileName = NULL;
}
#include <unistd.h>
SteppingAction::~SteppingAction()
{
  if(m_pFile!=NULL){
    fclose(m_pFile);
    unlink(m_szFileName);
  }
  if(m_szFileName) delete m_szFileName;
}

void SteppingAction::ForcedBeginEvent()
{
  m_dNi = 0;
}

#include <stdio.h>
#include <sys/mman.h>
#include <limits.h>
#include <fcntl.h>

G4double SteppingAction::GetSigmaFromFile()
{
  double Sigma=0;
  if(m_pFile!=NULL){
    fflush(m_pFile);
    fseek(m_pFile,0,SEEK_END);
    long size = ftell(m_pFile);
    fseek(m_pFile,0,SEEK_SET);
    struct __TmpFile{
      unsigned EvtID;
      double Weights;
    };
    if(size > UINT_MAX*sizeof(struct __TmpFile)){
      G4cout<<"Warning: file is too large. Some statistic will be truncated"<<G4endl;
      size = UINT_MAX*sizeof(struct __TmpFile);
    }
    unsigned fd = open(m_szFileName,0);
    struct __TmpFile* pTmpArr = (struct __TmpFile*)mmap(NULL,size,PROT_READ,MAP_PRIVATE,fd,0);
    if(pTmpArr==NULL){
      G4cout<<"WARNING: Error occur during maping the file into memory. No statistics in output"<<G4endl;
    }
    unsigned nEl = size/sizeof(struct __TmpFile);
    double Mean = 0;
    for(unsigned i=0;i<nEl;i++){
      Mean += pTmpArr[i].Weights/nEl;
    }
    Sigma = 0;
    for(i=0;i<nEl;i++){
      Sigma += pTmpArr[i].Weights*pTmpArr[i].Weights;
    }
    Sigma -= Mean*Mean;
    Sigma /= (nEl-1);
    munmap(pTmpArr,size);
    close(fd);
    fclose(m_pFile);
    unlink(m_szFileName);
    m_pFile = fopen(m_szFileName,"w");
  }
  return Sigma;
}

void SteppingAction::ForcedEndEvent()
{  
  if(m_pFile!=NULL){
    struct __TmpFile{
      unsigned EvtID;
      double Weights;
    }Output;
    Output.EvtID = m_pEvtAction->GetCounter();
    Output.Weights = m_dNi;
    fwrite(&Output,sizeof(Output),1,m_pFile);
  }
}
void SteppingAction::UseFile(char* szFileName)
{
  if(m_szFileName!=NULL) delete m_szFileName;
  if(m_pFile!=NULL) fclose(m_pFile);
  m_pFile = fopen(szFileName,"wb");
  if(m_pFile==NULL){
    m_szFileName=NULL;
    G4cout<<"Error opening file"<<szFileName<<G4endl;
    return;
  }
  m_szFileName = new char[strlen(szFileName)+1];
  strcpy(m_szFileName,szFileName);
}

void SteppingAction::UserSteppingAction(const G4Step* pStep)
{
  G4double dColimPos;
  G4Track* pTrack = pStep->GetTrack();
  G4ParticleDefinition* pParticle = pTrack->GetDefinition();
  G4DynamicParticle* pDynamic = (G4DynamicParticle*)pTrack->GetDynamicParticle();
  G4String szName = pParticle->GetParticleName();
  G4VPhysicalVolume* pVol = pTrack->GetVolume();
  G4String szVolName = (pVol!=NULL) ? pVol->GetName() : G4String("");
  G4VPhysicalVolume* pNextVol = pTrack->GetNextVolume();
  G4String szNextVolName = (pNextVol!=NULL) ? pNextVol->GetName() : G4String("");

  if(szName=="alpha") return;
  if(szName=="deuteron") return;
  if((szName!="proton")&&(szName!="neutron")){
    pTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  if(szName=="proton"){
    if(pDynamic->GetKineticEnergy()<7*keV) pTrack->SetTrackStatus(fStopAndKill);
    else if((szVolName == "Mother")&&(szNextVolName != "Target"))
      if(pTrack->GetParentID()==0) pTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  if((szName=="neutron")&&(pDynamic->GetKineticEnergy()<0.07*MeV)&&
     (szVolName!="Ghost")){
    pTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  else{
    dColimPos = m_pDetector->GetColimPos()/cm;
    dColimPos *= dColimPos;
    dColimPos *=5.94e-4;
    dColimPos = 1./(dColimPos*m_pRunAction->GetProtonsNumber());
    dColimPos *= pTrack->GetWeight();
    if(szVolName == "Detector"){
      if(pDynamic->GetMomentumDirection().z()>0){
	m_pAnalyzer->Fill(enHistoAtBeam,pDynamic->GetKineticEnergy()/eV,
			  dColimPos/m_pAnalyzer->GetEnergyRange(1,pDynamic->GetKineticEnergy()/eV)/eV/log(m_pAnalyzer->GetEnergy(1,pDynamic->GetKineticEnergy()/eV)*eV),pTrack->GetWeight());
	if(m_pEvtAction){
	  m_dNi += pTrack->GetWeight();
	}
      }
    }
    else if((szVolName=="Ghost")/*&&((szNextVolName!="Iron Shield A")&&(szNextVolName!="Iron Shield B"))*/){
      if((!m_pRunAction->GetStat())||(pTrack->GetParentID() == 0))
	if(pDynamic->GetMomentumDirection().z() > 0)
	  m_pAnalyzer->Fill(enHistoGhost,pDynamic->GetKineticEnergy()/MeV,pTrack->GetWeight()*1e+4/(5.94*sr*(pDynamic->GetKineticEnergy()/MeV)*m_pRunAction->GetProtonsNumber()),pTrack->GetWeight());
    }
    else if(szVolName=="Detector 20 cm"){
      if(pDynamic->GetMomentumDirection().z() > 0)
	m_pAnalyzer->Fill(enHisto20Cm,pDynamic->GetKineticEnergy()/eV,
			  dColimPos/m_pAnalyzer->GetEnergyRange(2,pDynamic->GetKineticEnergy()/eV)/eV/log(m_pAnalyzer->GetEnergy(2,pDynamic->GetKineticEnergy()/eV)*eV),pTrack->GetWeight());
    }
    else if(szVolName=="Detector 40 cm"){
      if(pDynamic->GetMomentumDirection().z()>0)
	m_pAnalyzer->Fill(enHisto40Cm,pDynamic->GetKineticEnergy()/eV,
			  dColimPos/m_pAnalyzer->GetEnergyRange(3,pDynamic->GetKineticEnergy()/eV)/eV/log(m_pAnalyzer->GetEnergy(3,pDynamic->GetKineticEnergy()/eV)*eV),pTrack->GetWeight());
    }
  }
}
