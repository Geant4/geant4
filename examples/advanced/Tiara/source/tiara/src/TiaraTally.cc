// $Id: TiaraTally.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "TiaraTally.hh"

TiaraTally::TiaraTally() :
  fBinEdges(),
  fLowestEdge(0),
  fMapHighEdgeToMeasure()
{}

TiaraTally::~TiaraTally()
{}

void TiaraTally::setBinEdges(const std::vector<double>  & binEdges){
  fBinEdges = binEdges;
  if (fBinEdges.size() < 2) {
    G4Exception("ERROR in  TiaraTally::setBinEdges: binEdges.size() < 2");
  }
  fLowestEdge = fBinEdges[0];
  std::vector<double>::iterator binIt = fBinEdges.begin();
  binIt++;
  for (;binIt != fBinEdges.end(); ++binIt) {
    std::vector<double>::iterator testIt = binIt;
    if (! (*binIt >  *(--testIt))) {
      G4Exception("ERROR in  TiaraTally::setBinEdges: binEdges not in ascending order");
    }

    std::pair<G4double, TiaraMeasure > myPair(0, TiaraMeasure());
    fMapHighEdgeToMeasure[*binIt] = myPair;
  }
}

void TiaraTally::fill(G4double x, G4double w){
  if (x>fLowestEdge) {
    for (TiaraTallyData::iterator dataIt = fMapHighEdgeToMeasure.begin();
	 dataIt != fMapHighEdgeToMeasure.end(); dataIt++) {
      if (x < dataIt->first) {
	dataIt->second.first += w;
	break;
      }
    }
  }
}

void TiaraTally::EndOfEventAction() {
  for (TiaraTallyData::iterator dataIt = fMapHighEdgeToMeasure.begin();
       dataIt != fMapHighEdgeToMeasure.end(); dataIt++) {
    dataIt->second.second.Xin(dataIt->second.first);
    dataIt->second.first = 0;
  }
}


TiaraMeasure TiaraTally::measure(G4int i){
  TiaraMeasure m;
  if (i<0 || i>(fMapHighEdgeToMeasure.size()-1)) {
    G4cout << "ERROR in TiaraTally::measure: argument: " << i 
	   << ", out of range!" <<  G4endl;
  } 
  else {
    G4int k(0);
    for (TiaraTallyData::iterator dataIt = fMapHighEdgeToMeasure.begin();
	 dataIt != fMapHighEdgeToMeasure.end(); dataIt++) {
      if (i==k++) {
	m = dataIt->second.second;
	break;
      }
    }
  }
  return m;
}

std::vector<double>  TiaraTally::binEdges(){
  return fBinEdges;
}

G4int TiaraTally::size() {
  return fMapHighEdgeToMeasure.size();
}







