#ifndef MLDoseAnalyser_h
#define MLDoseAnalyser_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#include "MLAnalysisManager.hh"
#include "MLGeometryConstruction.hh"

#ifndef USEHBOOK
  #include "CSVofstream.hh"
  #include "RPTofstream.hh"
#endif
////////////////////////////////////////////////////////////////////////////////
//
class MLDoseAnalyser
{
public:
  MLDoseAnalyser (MLAnalysisManager*, MLGeometryConstruction*);
  ~MLDoseAnalyser ();

  void SetDoseUnit (G4String s) {doseUnit = s;};

  void AddDoseLayer (G4int);
  void DeleteDoseLayer (G4int);
  void ListDoseLayer ();
  void DeleteLayer (G4int);

  G4int GetNBlocks ();
  G4int GetAddedNBlocks (G4int);

  void BeginOfRunAction (G4double);
  void AddToDoseDetector (G4int, G4double);

#ifndef USEHBOOK
  void NormaliseOutput (G4double);
#endif
  G4int GetDLayerIdx (G4int i) {return doseLayers[i];};
  G4int GetNbOfDLayers () {return G4int(doseLayers.size());};
  G4double GetDose (G4int i) {return dose[i];};
  G4double GetDoseErr (G4int i) {return doseErr[i];};
  G4String GetDoseUnit () {return doseUnit;};

private:
  std::vector<G4double>  ADose;
  std::vector<G4double>  ADoseSqr;
  std::vector<G4double>  dose;
  std::vector<G4double>  doseErr;
  std::vector<G4int>     doseLayers;
  G4String                 doseUnit;

  std::vector<G4double>  factor;
  std::vector<G4double>  factorSqr;
  MLAnalysisManager       *analysisManager;
  MLGeometryConstruction  *geometry;
  MLMaterial*              materialsManager;
  G4int                    Nb;

#ifndef USEHBOOK
  friend CSVofstream & operator << (CSVofstream &, MLDoseAnalyser &);
  friend RPTofstream & operator << (RPTofstream &, MLDoseAnalyser &);
#endif
};
////////////////////////////////////////////////////////////////////////////////
#endif
