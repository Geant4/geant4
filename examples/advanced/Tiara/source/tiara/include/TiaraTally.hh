#ifndef TiaraTally_hh
#define TiaraTally_hh TiaraTally_hh
#include "TiaraMeasure.hh"
#include <vector>
#include "g4std/map"


typedef G4std::map<G4double, G4std::pair<G4double, TiaraMeasure > > TiaraTallyData;

class TiaraTally{
public:
  TiaraTally();
  ~TiaraTally();
	     
  void setBinEdges(const std::vector<double>  & binEdges);
  void fill(G4double x, G4double w);
  void EndOfEventAction();
  TiaraMeasure measure(G4int i);
  std::vector<double>  binEdges();
  G4int size();

private:
  std::vector<double>  fBinEdges;
  G4double fLowestEdge;
  TiaraTallyData fMapHighEdgeToMeasure;
  
  

};

#endif
