#ifndef MLPHSAnalyser_h
#define MLPHSAnalyser_h 1
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

class MLAnalysisManager;
class MLGeometryConstruction;
////////////////////////////////////////////////////////////////////////////////
//
class MLPHSAnalyser
{
public:
  MLPHSAnalyser (MLAnalysisManager*, MLGeometryConstruction*);
  ~MLPHSAnalyser ();

  void SetPType (G4String);
  void SetPHSMax (G4double);
  void SetPHSMin (G4double);
  void SetPNBin (G4int);
  G4int GetPNBin () {return pNBin;}
  void AddPEdg (G4double);
  void DeletePEdg (G4double);
  void ClearPEdg () {pEdg.clear();};
  void ListPEdg ();
  void DefaultPEdg ();

  G4int GetNBlocks ();
  G4int GetAddedNBlocks (G4int);

  void BeginOfRunAction (G4double);
  void FillHbook (G4int, G4double, G4double, G4double);
#ifndef USEHBOOK
  void NormaliseOutput (G4double);
#endif

#ifdef USEHBOOK
  HbookHistogram* GetHisto (G4int i) {return Histo[i];};
#else
  MLHisto1D* GetHisto (G4int i) {return Histo[i];};
#endif

private:
  ScaleType                   pType;
  G4double                    maxphs;
  G4double                    minphs;
  G4int                       pNBin;
  std::vector<G4double>     pEdg;

#ifdef USEHBOOK
  std::vector<HbookHistogram*>  Histo;
#else
  std::vector<MLHisto1D*>   Histo;
  CSVofstream                 CSVfile;
  RPTofstream                 RPTfile;
#endif

  G4double                    factor;
  G4int                       Nb;
  MLAnalysisManager          *analysisManager;
  MLGeometryConstruction     *geometry;

#ifndef USEHBOOK
  friend CSVofstream & operator << (CSVofstream &, MLPHSAnalyser &);
  friend RPTofstream & operator << (RPTofstream &, MLPHSAnalyser &);
#endif
};
////////////////////////////////////////////////////////////////////////////////
#endif
