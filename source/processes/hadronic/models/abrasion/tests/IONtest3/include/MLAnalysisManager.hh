#ifndef MLAnalysisManager_h
#define MLAnalysisManager_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#ifdef USEHBOOK
  #include "HbookHistogram.hh"
#else
  #include "CSVofstream.hh"
  #include "RPTofstream.hh"
  #include "MLHisto1D.hh"
#endif

class MLGeometryConstruction;
class MLAnalysisMessenger;
class MLFluenceAnalyser;
class MLDoseAnalyser;
class MLNIELAnalyser;
class MLPHSAnalyser;
////////////////////////////////////////////////////////////////////////////////
//
class MLAnalysisManager
{
private:
  MLAnalysisManager (int = 0, char** = 0);

public:
  virtual ~MLAnalysisManager ();
  static MLAnalysisManager* getInstance (int = 0, char** = 0);

  MLFluenceAnalyser* GetFluenceAnalyser () {return fluenceAnalyser;};
  MLDoseAnalyser*    GetDoseAnalyser ()    {return doseAnalyser;};
  MLNIELAnalyser*    GetNIELAnalyser ()    {return NIELAnalyser;};
  MLPHSAnalyser*     GetPHSAnalyser ()     {return PHSAnalyser;};

  void BeginOfRunAction ();
  void EndOfRunAction (G4double);
  void EndOfEventAction (G4double);

  void SetFilename (G4String file) {filenamePrefix = file;};
  void SetNormalFactor (G4double d) {nFactor = d;};
  const G4double GetNormalFactor () {return nFactor;};

#ifndef USEHBOOK
private:
  void NormaliseOutput (G4double);
  void PrintOutput (G4double);
#endif


private:
  static MLAnalysisManager  *instance;

  G4double                   nFactor;
  MLGeometryConstruction    *geometry;
  MLAnalysisMessenger       *analysisMessenger;

  G4String                   filenamePrefix;
  MLFluenceAnalyser         *fluenceAnalyser;
  MLDoseAnalyser            *doseAnalyser;
  MLNIELAnalyser            *NIELAnalyser;
  MLPHSAnalyser             *PHSAnalyser;

#ifdef USEHBOOK
  std::vector<HbookHistogram*>  Histo;
#else
  std::vector<MLHisto1D*>  Histo;
  CSVofstream                CSVFile;
  RPTofstream                RPTFile;
#endif

};
////////////////////////////////////////////////////////////////////////////////
#endif
