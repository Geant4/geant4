//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TiaraTally.cc,v 1.5 2006/06/29 15:45:37 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
  if (i<0 || i>(G4int)(fMapHighEdgeToMeasure.size()-1)) {
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







