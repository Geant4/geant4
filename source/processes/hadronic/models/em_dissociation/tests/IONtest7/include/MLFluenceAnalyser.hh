#ifndef MLFluenceAnalyser_h
#define MLFluenceAnalyser_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#include "MLScaleType.hh"

#ifdef USEHBOOK
  #include "HbookHistogram.hh"
#else
  #include "CSVofstream.hh"
  #include "RPTofstream.hh"
  #include "MLHisto1D.hh"
#endif

class MLGeometryConstruction;
class MLAnalysisManager;
////////////////////////////////////////////////////////////////////////////////
//
class MLFluenceAnalyser
{
public:
  MLFluenceAnalyser (MLAnalysisManager*, MLGeometryConstruction*);
  ~MLFluenceAnalyser ();

  void SetFluenceType (G4String);
  G4String GetFluenceType ();
  G4bool IsDividedByCosTheta () {return divideByCosT;}
  void SetFluenceUnit (G4String s) {fluenceUnit = s;};
  G4String GetFluenceUnit () {return fluenceUnit;};

  void SetEType (G4String);
  void SetEngMax (G4double max);
  void SetEngMin (G4double min);
  void SetENBin (G4int);
  G4int GetENBin () {return eNBin;}
  void AddEEdg (G4double);
  void DeleteEEdg (G4double);
  void ClearEEdg () {eEdg.clear();};
  void ListEEdg ();
  void DefaultEEdg ();

  void SetAType (G4String);
  void SetAngMax (G4double max);
  void SetAngMin (G4double min);
  void SetANBin (G4int);
  G4int GetANBin () {return aNBin;}
  void AddAEdg (G4double);
  void DeleteAEdg (G4double);
  void ClearAEdg () {aEdg.clear();};
  void ListAEdg ();
  void DefaultAEdg ();
  G4double GetAEdg (G4int i) {return aEdg[i];};

  void AddSPart (G4String);
  void DeleteSPart (G4String);
  G4bool IsParticleOutput (G4String part)
    {return binary_search( sPart.begin(), sPart.end(), part );}
  void ClearSPart () {sPart.clear();};
  void ListSPart ();

  void BeginOfRunAction (G4double);
  G4int GetHistIdx (G4int i, G4int j, G4int k);
  void TallyFluenceEvent (G4int, G4int, G4double, G4double, G4double);

  G4int GetNBlocks ();
  G4int GetAddedNBlocks (G4int);
#ifndef USEHBOOK
  void NormaliseOutput (G4double);
#endif
  void EndOfRunAction ();

#ifdef USEHBOOK
  HbookHistogram* GetHisto (G4int i) {return Histo[i];};
#else
  MLHisto1D* GetHisto (G4int i) {return Histo[i];};
#endif

private:
  void FillHbook (G4int, G4double, G4double, G4double);

private:
  std::vector<G4double>  factor;
  MLAnalysisManager       *analysisManager;
  MLGeometryConstruction  *geometry;
  G4bool                   divideByCosT;
  G4String                 fluenceUnit;

  ScaleType                eType;
  G4double                 maxeng;
  G4double                 mineng;
  G4int                    eNBin;
  std::vector<G4double>  eEdg;
  
  ScaleType                aType;
  G4double                 maxang;
  G4double                 minang;
  G4int                    aNBin;
  std::vector<G4double>  aEdg;

  std::vector<G4String>  sPart;

#ifdef USEHBOOK
  std::vector<HbookHistogram*> Histo;
#else
  std::vector<MLHisto1D*> Histo;
#endif

  G4int                    Nb;

#ifndef USEHBOOK
  friend CSVofstream & operator << (CSVofstream &, MLFluenceAnalyser &);
  friend RPTofstream & operator << (RPTofstream &, MLFluenceAnalyser &);
#endif

};
////////////////////////////////////////////////////////////////////////////////
#endif
