#ifdef G4ANALYSIS_USE
#ifndef fluoTestAnalysisManager_h
#define fluoTestAnalysisManager_h 1

#include "G4VAnalysisManager.hh"

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class fluoTestAnalysisMessenger;
class fluoTestDetectorConstruction;
class IHistogramFactory;
class IHistogram1D;
class IPlotter;
class IVectorFactory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestAnalysisManager: public G4VAnalysisManager
{
public:
  fluoTestAnalysisManager(fluoTestDetectorConstruction*);
  virtual ~fluoTestAnalysisManager();
  
public:
  G4bool RegisterAnalysisSystem(G4VAnalysisSystem*);
  IHistogramFactory* GetHistogramFactory(const G4String&);
 
  void Store(IHistogram* = 0, const G4String& = "");
  void Plot(IHistogram* = 0);
  void InsertKEnergy(double gKe);
  void InsertIoniEnergy(double Ioe);
  void InsertPhotoEnergy(double Phe);
  void InsertBremEnergy(double Bre);
  void InsertComptEnergy(double Coe);
  void InsertConvEnergy(double Cve);
  void InsertRaylEnergy(double Rae);
  void InsertOutEnergy(double Goe);
 
  void InserteKEnergy(double eKe);
  void InserteIoniEnergy(double eIoe);
  void InsertePhotoEnergy(double ePhe);
  void InserteBremEnergy(double eBre);
  void InserteComptEnergy(double eCoe);
  void InserteConvEnergy(double eCve);
  void InserteRaylEnergy(double eRae);
  void InserteOutEnergy(double Eoe);

 void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void SetHisto1DDraw(G4String str) {histo1DDraw = str;};
  void SetHisto1DSave(G4String str) {histo1DSave = str;};
 
private:
  G4VAnalysisSystem* analysisSystem;
  IPlotter* pl;
  IVectorFactory* fVectorFactory;
  IHistogramFactory* histoFactory;

  IHistogram1D*  histoGammaKenergy;
  IHistogram1D* histoIoniEnergy;
  IHistogram1D*  histoPhotoEnergy;
  IHistogram1D* histoBremEnergy;
  IHistogram1D*  histoComptEnergy;
  IHistogram1D*  histoConvEnergy;
  IHistogram1D*  histoRaylEnergy;
  IHistogram1D*  histoGammaOutEnergy;
  IHistogram1D*  histoElecKenergy;
  IHistogram1D* histoeIoniEnergy;
  IHistogram1D*  histoePhotoEnergy;
  IHistogram1D* histoeBremEnergy;
  IHistogram1D*  histoeComptEnergy;
  IHistogram1D*  histoeConvEnergy;
  IHistogram1D*  histoeRaylEnergy;
  IHistogram1D*  histoElecOutEnergy;
  fluoTestDetectorConstruction*    Detector;

  G4String histo1DDraw;
  G4String histo1DSave;
 

  fluoTestAnalysisMessenger* analysisMessenger;
};


#endif
#endif



