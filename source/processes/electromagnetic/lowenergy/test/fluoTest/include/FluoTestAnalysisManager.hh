#ifdef G4ANALYSIS_USE
#ifndef FluoTestAnalysisManager_h
#define FluoTestAnalysisManager_h 1

#include "G4VAnalysisManager.hh"

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class FluoTestAnalysisMessenger;
class FluoTestDetectorConstruction;
class IHistogramFactory;
class IHistogram1D;
class IPlotter;
class IVectorFactory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestAnalysisManager: public G4VAnalysisManager
{
public:
  FluoTestAnalysisManager(FluoTestDetectorConstruction*);
  virtual ~FluoTestAnalysisManager();
  
public:
  G4bool RegisterAnalysisSystem(G4VAnalysisSystem*);
  IHistogramFactory* GetHistogramFactory(const G4String&);
 
  void Store(IHistogram* = 0, const G4String& = "");
  void Plot(IHistogram* = 0);
  void InsGamBornSample(double gBs);
  void InsEleBornSample(double eBs);
  // void InsGamLSBackw(double gSb);
  void InsGamLeavSam(double gLs);
  void InsEleLeavSam(double eLs);
  // void InsGamDetPre(double gDpr);
  // void InsGamDetPost(double gDps);
  //  void InsOtherPart(double oP);
  void InsDetETot(double dEt);
  void InsGamDet(double gD);
  
  void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);
  
  void SetHisto1DDraw(G4String str) {histo1DDraw = str;};
  void SetHisto1DSave(G4String str) {histo1DSave = str;};
 
private: 
  
  FluoTestDetectorConstruction*    Detector;

  //IHistogram1D*  histoGamDetPre;
  // IHistogram1D* histoGamDetPost;
  IHistogram1D*  histoGamLeavSam;
  IHistogram1D*  histoEleLeavSam;
  //IHistogram1D* histoGamLSBack;
  IHistogram1D*  histoGamBornSam;
  IHistogram1D*  histoEleBornSam;
  //IHistogram1D*  histoOtherPartDet;
  //IHistogram1D*  histoDetETot;
  //IHistogram1D* histoGamDet;
   IHistogramFactory* histoFactory;
  IPlotter* pl;
  G4VAnalysisSystem* analysisSystem;
  IVectorFactory* fVectorFactory;

  G4String histo1DDraw;
  G4String histo1DSave;
 

  FluoTestAnalysisMessenger* analysisMessenger;
};


#endif
#endif



