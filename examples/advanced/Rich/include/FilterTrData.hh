// Rich advanced example for Geant4
// FilterTrData.chh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef FilterTrData_h
#define FilterTrData_h 1
#include "G4ios.hh" 
#include "globals.hh"
#include "g4std/vector"
#include "FilterTypeSpec.hh"

class FilterTrData {

public: 
  FilterTrData(G4int , G4String );
  virtual ~FilterTrData();

  G4double GetCurrentFilTrans(G4double);
  G4double GetTotFilTrans(G4double);

  void ReadFilterTrans();
  G4int GetNumberofBins() {return NumberOfTrBins; }
  G4double GetCurTransWL(G4int BinNumw ) {return TransWaveL[BinNumw]; }
  G4double GetCurTransValue(G4int BinNumv ) {return TransValue[BinNumv];}
  G4double GetCurTransTotValue(G4int BinNumv ) {return TransTotValue[BinNumv];}
  
  G4double GetCurFilterThickness() {return FilterThickness; }
  vector<G4double> GetTransWL() { return TransWaveL ; }
  vector<G4double> GetTransValue() {return TransValue ; }
  vector<G4double> GetTransTotValue() {return TransTotValue ; }
  FilterType GetFilterTypeIndex() {return FilterTypeIndex ; }

private:

  FilterType FilterTypeIndex;
  G4int FilterTypeNumber;
  G4int NumberOfTrBins;
  G4double FilterThickness;
  G4double FilterRefIndexNominal;
  vector<G4double>TransWaveL ;
  vector<G4double>TransValue ;
  vector<G4double>TransTotValue;
  G4double CurNeighbourRefIndexNominal;
  G4String FilterTransDataFileName;

};

#endif 
