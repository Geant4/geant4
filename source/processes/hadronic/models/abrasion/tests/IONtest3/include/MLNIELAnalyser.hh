#ifndef MLNIELAnalyser_h
#define MLNIELAnalyser_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#include "MLNIELType.hh"

#ifdef USEHBOOK
  #include "HbookHistogram.hh"
#else
  #include "CSVofstream.hh"
  #include "RPTofstream.hh"
  #include "MLHisto1D.hh"
#endif

#include "MLNIELFunction.hh"

class MLGeometryConstruction;
class MLFluenceAnalyser;
class MLAnalysisManager;
////////////////////////////////////////////////////////////////////////////////
//
class MLNIELAnalyser
{
public:
  MLNIELAnalyser (MLAnalysisManager*, MLGeometryConstruction*);
  ~MLNIELAnalyser ();

  void SelectNielFunction (G4String);
  NIELType GetNielType () {return nielType;};
  void SetDoseUnit (G4String s) {doseUnit = s;};

  void AddNielLayer (G4int);
  void DeleteNielLayer (G4int);
  void ListNielLayer ();
  void DeleteLayer (G4int);
  G4int GetNbOfNielLayers () {return nielLayers.size();};
  G4int GetNLayerIdx (G4int i) {return nielLayers[i];}

  void BeginOfRunAction (G4double);
  G4int GetNBlocks ();
  G4int GetAddedNBlocks (G4int);
  void BeginOfEventAction ();
  void EndOfEventAction ();
  void AddToNIELDetector (G4int, G4int, G4double, G4double, G4double);

#ifndef USEHBOOK
  void NormaliseOutput (G4double);
#endif

  G4double GetTotal       (G4int i) {return totalniel[i];}
  G4double GetTotalErr    (G4int i) {return totalErr[i];}
  G4double GetProton      (G4int i) {return protonniel[i];}
  G4double GetProtonErr   (G4int i) {return protonErr[i];}
  G4double GetNeutron     (G4int i) {return neutronniel[i];}
  G4double GetNeutronErr  (G4int i) {return neutronErr[i];}
  G4double GetElectron    (G4int i) {return electronniel[i];}
  G4double GetElectronErr (G4int i) {return electronErr[i];}
  G4double GetPion        (G4int i) {return pionniel[i];}
  G4double GetPionErr     (G4int i) {return pionErr[i];}
  G4String GetDoseUnit () {return doseUnit;};

private:

  std::vector<G4double> totalniel;
  std::vector<G4double> totalErr;
  std::vector<G4double> protonniel;
  std::vector<G4double> protonErr;
  std::vector<G4double> neutronniel;
  std::vector<G4double> neutronErr;
  std::vector<G4double> electronniel;
  std::vector<G4double> electronErr;
  std::vector<G4double> pionniel;
  std::vector<G4double> pionErr;

  std::vector<G4double> ATotal;
  std::vector<G4double> ATotalSqr;
  std::vector<G4double> AProton;
  std::vector<G4double> AProtonSqr;
  std::vector<G4double> ANeutron;
  std::vector<G4double> ANeutronSqr;
  std::vector<G4double> AElectron;
  std::vector<G4double> AElectronSqr;
  std::vector<G4double> APion;
  std::vector<G4double> APionSqr;

  std::vector<G4double> EvtProtonniel;
  std::vector<G4double> EvtNeutronniel;
  std::vector<G4double> EvtElectronniel;
  std::vector<G4double> EvtPionniel;

  std::vector<G4int>    nielLayers;

#ifdef USEHBOOK
  std::vector<HbookHistogram*> Histo;
#else
  std::vector<MLHisto1D*> Histo;
#endif

  std::vector<G4double>  factor;
  std::vector<G4double>  factorSqr;
  NIELType                 nielType;
  MLAnalysisManager       *analysisManager;
  MLFluenceAnalyser       *fluenceAnalyser;
  MLGeometryConstruction  *geometry;
  MLNIELFunction          *nielFunc;
  G4String                 doseUnit;
  G4int                    Nb;

#ifndef USEHBOOK
  friend CSVofstream & operator << (CSVofstream &, MLNIELAnalyser &);
  friend RPTofstream & operator << (RPTofstream &, MLNIELAnalyser &);
#endif
};
////////////////////////////////////////////////////////////////////////////////
#endif
