// $Id: TiaraTally.hh,v 1.3 2003-06-18 16:40:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraTally
//

#ifndef TiaraTally_hh
#define TiaraTally_hh TiaraTally_hh
#include "TiaraMeasure.hh"
#include <vector>
#include <map>


typedef std::map<G4double, std::pair<G4double, TiaraMeasure > > TiaraTallyData;

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
