//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: TiaraTally.hh,v 1.4 2003/06/25 09:12:52 gunter Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
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
