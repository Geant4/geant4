#ifdef  G4ANALYSIS_USE

#include <stdlib.h>
#include "g4std/fstream"
#include "FluoTestAnalysisManager.hh"
#include "G4VAnalysisSystem.hh"
#include "FluoTestDetectorConstruction.hh"
#include "FluoTestAnalysisMessenger.hh"
#include <IHistogramFactory.h>
#include <IHistogram1D.h>
#include "G4ios.hh"
#include <IPlotter.h>
#include <IVector.h>
#include <IVectorFactory.h>
#include "g4rw/tvordvec.h"
#include "G4LizardSystem.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisManager::FluoTestAnalysisManager(FluoTestDetectorConstruction* FluoTestDC): 
  Detector(FluoTestDC),
 
  //histoGamDet(0),
  //histoDetETot(0),
  //histoGamDetPre(0),
  //histoGamDetPost(0),
  histoGamLeavSam(0),
  histoEleLeavSam(0),
  //  histoGamLSBack(0),
  histoGamBornSam(0),
  histoEleBornSam(0),
  // histoOtherPartDet(0),
  
  histoFactory(0), pl(0),
  histo1DDraw("disenable"),histo1DSave("enable")
  
{
  // Define the messenger and the analysis system
  analysisMessenger = new FluoTestAnalysisMessenger(this);
  analysisSystem = new G4LizardSystem;

  histoFactory = analysisSystem->GetHistogramFactory();  

  //   The following lines set the plotter and the vectorfactory that
  //   are needed in this example for a multiple histograms
  //   visualization.
 
 fVectorFactory = createIVectorFactory();          
  pl = createIPlotter();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisManager::~FluoTestAnalysisManager() {

  //delete histoGamDet;
  //delete histoDetETot;
  // delete histoGamDetPre;
  //delete histoGamDetPost;
  delete histoGamLeavSam;
  delete histoEleLeavSam;
  //    delete histoGamLSBack; 
  delete histoGamBornSam;
  delete histoEleBornSam;
  // delete histoOtdelete 
  
  delete analysisSystem;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FluoTestAnalysisManager::RegisterAnalysisSystem(G4VAnalysisSystem*)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

IHistogramFactory* FluoTestAnalysisManager::GetHistogramFactory(const G4String& aSystem)
{
  return histoFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestAnalysisManager::Store(IHistogram* histo, const G4String& ID)
{
  analysisSystem->Store(histo, ID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void FluoTestAnalysisManager::Plot(IHistogram* histo = 0)
{ 
  // In a normal case we use the following line to use the standard plot
  // analysisSystem->Plot(histo);
  
  //IVector* vgD = 0;
  //IVector* vdEt = 0;
  IVector* vgBs = 0;
  IVector* veBs = 0;
  //  IVector* vgSb = 0;
  IVector* vgLs = 0;
  IVector* veLs = 0;
  IVector* vgDpr = 0;
  // IVector* vgDps = 0;
  // IVector* voP = 0;
  
  // We fill them with the histograms
  
  if(histo1DDraw == "enable")
    {
     
      // vgD = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamDet));
      //vdEt = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoDetETot));
      
      vgBs = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamBornSam));
      veBs = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoEleBornSam));
      //      vgSb = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamLSBack));    
       vgLs = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamLeavSam));
       veLs = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoEleLeavSam));  
       // vgDpr = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamDetPre));
       
       
      // vgDps = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamDetPost));
       // voP = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoOtherPartDet));
       // pl->plot(vgD);
       //  pl->plot(vdEt);
       pl->plot(vgBs);
       pl->plot(veBs);
       //   pl->plot(vgSb); 
       pl->plot(vgLs);
       pl->plot(veLs);
       // pl->plot(vgDpr); 
       // pl->plot(vgDps);
       // pl->plot(voP); 
       
       pl->refresh(); 
    }
  //delete vgD;
  //delete vdEt;
 delete vgBs;
 delete veBs;
 // delete vgSb;
 delete vgLs;
 delete veLs;
 delete vgDpr;
 // delete vgDps;
 //delete voP;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/*
  void FluoTestAnalysisManager::InsGamDet(double gD)
  {
  histoGamDet->fill(gD);
  }
 
  void FluoTestAnalysisManager::InsDetETot(double dEt)
  {
  histoDetETot->fill(dEt);
  }
*/
void FluoTestAnalysisManager::InsGamBornSample(double gBs)
{
  histoGamBornSam->fill(gBs);
}
void FluoTestAnalysisManager::InsEleBornSample(double eBs)
{
  histoEleBornSam->fill(eBs);
}
/*
  void FluoTestAnalysisManager::InsGamLSBackw(double gSb)
{
histoGamLSBack->fill(gSb);
}
*/
void FluoTestAnalysisManager::InsGamLeavSam(double gLs)
{
  histoGamLeavSam->fill(gLs);
}
void FluoTestAnalysisManager::InsEleLeavSam(double eLs)
{
  histoEleLeavSam->fill(eLs);
}
/*
  void FluoTestAnalysisManager::InsGamDetPre(double gDpr)
  {
   histoGamDetPre->fill(gDpr);
   }*/
/*

  void FluoTestAnalysisManager::InsGamDetPost(double gDps)
  {
  histoGamDetPost->fill(gDps);
  }
  
void FluoTestAnalysisManager::InsOtherPart(double oP)
{
histoOtherPartDet->fill(oP);
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This member reset the histograms and it is called at the begin
// of each run; here we put the inizialization so that the histograms have 
// always the right dimensions depending from the detector geometry

void FluoTestAnalysisManager::BeginOfRun() 
{ 
  if (histoFactory) 
    {
      /*   histoFactory->destroy(histoGamDet);
	   histoGamDet = histoFactory->create1D("Energy released by the gammas(keV)", 500, 0,7);
	   
	   histoFactory->destroy(histoDetETot);
      histoDetETot = histoFactory->create1D("Total energy in the detector(keV)", 500, 0,7);
      histoFactory->destroy(histoGamDetPre);
      histoGamDetPre = histoFactory->create1D("Gammas at detector pre step point(keV)", 500, 0,5);*/
      /*
	histoFactory->destroy(histoGamDetPost);
	histoGamDetPost = histoFactory->create1D("Gammas at detector post step point(keV)", 500, 0,1000);
      */
      histoFactory->destroy(histoGamLeavSam);
      histoGamLeavSam = histoFactory->create1D("Gammas leaving sample(keV)", 500, 0,6);
      histoFactory->destroy(histoEleLeavSam);
      histoEleLeavSam = histoFactory->create1D("Electrons leaving sample(keV)", 500, 0,6);
      
      /*  
	  histoFactory->destroy(histoGamLSBack);
	  histoGamLSBack = histoFactory->create1D("Gammas leaving sample backward(keV)", 500, 0,1000);*/
      
      histoFactory->destroy(histoGamBornSam);
      histoGamBornSam = histoFactory->create1D("Gammas generated in sample (keV)", 500, 0,6);
      
      histoFactory->destroy(histoEleBornSam);
      histoEleBornSam = histoFactory->create1D("Electrons generated in sample (keV)", 500, 0,6);
	/*
	  histoFactory->destroy(histoOtherPartDet);
	  histoOtherPartDet = histoFactory->create1D("Other particles reaching detector(keV)", 500, 0,1000);
      */        
    }
  /*
    if(histoDetETot)
    histoDetETot->reset();
    
    if(histoGamDetPre)
    histoGamDetPre->reset();
  */  
  /*
  if( histoGamDetPost)
  histoGamDetPost ->reset();
  */ 
  if(histoGamLeavSam)
    histoGamLeavSam ->reset();
  if(histoEleLeavSam)
    histoEleLeavSam ->reset();
  // if(histoGamLSBack)
  //    histoGamLSBack ->reset();
  
  if(histoGamBornSam)
    histoGamBornSam ->reset();
  if(histoEleBornSam)
    histoGamBornSam ->reset();
  /*   if (histoOtherPartDet)
       histoOtherPartDet ->reset();
  */
  
  pl->zone(1,1);
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
//   This member is called at the end of each run 

void FluoTestAnalysisManager::EndOfRun(G4int n) 
{
  // This variable contains the names of the PS files
  char name[15];
  // We define some vectors
  //  IVector* vgD = 0;
  // IVector* vdEt = 0;
  IVector* vgBs   = 0;
  IVector* veBs   = 0;
  //  IVector* vgSb   = 0;
  IVector* vgLs   = 0;
  IVector* veLs   = 0;
  //IVector* vgDpr   = 0;
  //IVector* vgDps   = 0;
  //IVector* voP   = 0;
  
  // Temporary we set one single zone for the plotter
  pl->zone(1,1);
  
  // We now print the histograms, each one in a separate file 
  
  if(histo1DSave == "enable")
    {
      //    vgD = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamDet));
      // vgD->toAscii("gammaDet.dat");
      // vdEt=fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoDetETot));
      vgBs   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamBornSam));
      veBs   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoEleBornSam)); 
      //  vgSb   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamLSBack));
      vgLs   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamLeavSam));
      veLs   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoEleLeavSam));
      //vgDpr   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamDetPre));
      //vgDps  = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoGamDetPost));
      //voP   = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(histoOtherPartDet));
      /*   
      sprintf(name,"GamDet.ps", n);
      pl->plot(vgD);
      pl->psPrint(name);
      
      sprintf(name,"TotEnDet.ps", n);
      pl->plot(vdEt);
      pl->psPrint(name);*/

      //sprintf(name,"gamDetPre.ps", n);
      //pl->plot(vgDpr);
      //pl->psPrint(name);
      /*
	sprintf(name,"gamDetPost.ps", n);
	pl->plot(vgDps);
	pl->psPrint(name);
      */
      sprintf(name,"gamLeavSam.ps", n);
      pl->plot(vgLs);
      pl->psPrint(name);
      sprintf(name,"eleLeavSam.ps", n);
      pl->plot(veLs);
      pl->psPrint(name);
      /*  
	  sprintf(name,"gamLeavSBack.ps", n);
	  pl->plot(vgSb);
	  pl->psPrint(name);
      */
      // vgBs->toAscii("gamBornSam.dat");
      sprintf(name,"gamBornSam.ps", n);
      pl->plot(vgBs);
      pl->psPrint(name);
      
      sprintf(name,"eleBornSam.ps", n);
      pl->plot(veBs);
      pl->psPrint(name);
      
      /*
	sprintf(name,"othPartDet.ps", n);
	pl->plot(voP);
	pl->psPrint(name);
      */
      
    }
  //delete vgD;
  //delete vdEt;
  //delete vgDpr;
  //delete vgDps;
  delete vgLs;
  delete veLs;
  // delete vgSb;
  delete vgBs;
  delete veBs;
  
  //delete voP;
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//  This member is called at the end of every event 
void FluoTestAnalysisManager::EndOfEvent(G4int flag) 
{
  // The histograms are updated only if there is some
  // hits in the event
  if(flag) Plot();
}

#endif







