#ifdef G4ANALYSIS_USE
#ifndef FluoTestAnalysisManager_h
#define FluoTestAnalysisManager_h 1

//#include "G4VAnalysisManager.hh"

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class FluoTestAnalysisMessenger;
class FluoTestDetectorConstruction;
//class IHistogramFactory;
class IHistoManager;
class IHistogram1D;
//class IPlotter;
namespace Lizard {
class NTupleFactory;
class NTuple;
};
//class IVectorFactory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestAnalysisManager
//: public G4VAnalysisManager
{
public:
  FluoTestAnalysisManager(FluoTestDetectorConstruction*);
  virtual ~FluoTestAnalysisManager();
  
public:
  
  void InsGamBornSample(double gBs);
  void InsEleBornSample(double eBs);
   void InsGamLS(double gS);
 void InsGamLSP(double gSP);
  void InsGamLeavSam(double gLs);
  void InsEleLeavSam(double eLs);
   void InsGamDetPre(double gDpr);
  // void InsGamDetPost(double gDps);
  //  void InsOtherPart(double oP);
  void InsDetETot(double dEt);
  void InsGamDet(double gD);
  
  void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);
  
 
private:
  
  IHistoManager* histoManager;
  Lizard::NTupleFactory* factory;

  FluoTestDetectorConstruction*    Detector;
  Lizard::NTuple* ntuple;
  IHistogram1D* histoGamDet;
  IHistogram1D*  histoGamDetPre;
  // IHistogram1D* histoGamDetPost;
  IHistogram1D*  histoGamLeavSam;
  IHistogram1D*  histoEleLeavSam;
  IHistogram1D* histoGamLS;
 IHistogram1D* histoGamLSP;
  IHistogram1D*  histoGamBornSam;
  IHistogram1D*  histoEleBornSam;
  //IHistogram1D*  histoOtherPartDet;
  //IHistogram1D*  histoDetETot;
 
  FluoTestAnalysisMessenger* analysisMessenger;
};


#endif
#endif



