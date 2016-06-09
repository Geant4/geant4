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
//
// $Id: G4WeightWindowAlgorithm.cc,v 1.7 2003/08/19 15:16:56 dressel Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WeightWindowAlgorithm.cc
//
// ----------------------------------------------------------------------
#include "G4WeightWindowAlgorithm.hh"
#include "Randomize.hh"


G4WeightWindowAlgorithm::
G4WeightWindowAlgorithm(G4double upperLimitFaktor,
			G4double survivalFaktor,
			G4int maxNumberOfSplits) :
  fUpperLimitFaktor(upperLimitFaktor),
  fSurvivalFaktor(survivalFaktor),
  fMaxNumberOfSplits(maxNumberOfSplits)
{}

G4WeightWindowAlgorithm::~G4WeightWindowAlgorithm()
{}



G4Nsplit_Weight 
G4WeightWindowAlgorithm::Calculate(G4double init_w,
				   G4double lowerWeightBound) const {

  G4double survivalWeight = lowerWeightBound * fSurvivalFaktor;
  G4double upperWeight = lowerWeightBound * fUpperLimitFaktor;

  // initialize return value for case weight in window
  G4Nsplit_Weight nw;
  nw.fN = 1;
  nw.fW = init_w;

  if (init_w > upperWeight) {
    // splitting

    G4double wi_ws = init_w/survivalWeight;
    G4int int_wi_ws = static_cast<int>(wi_ws);

    // values in case integer mode or in csae of double
    // mode and the lower number of splits will be diced
    nw.fN = int_wi_ws;
    nw.fW = survivalWeight;	

    if (wi_ws <= fMaxNumberOfSplits) {
      if (wi_ws > int_wi_ws) {
	// double mode
	G4double p2 =  wi_ws - int_wi_ws;
	G4double r = G4UniformRand();
	if (r<p2) {
	  nw.fN = int_wi_ws + 1;
	}
      }
    }
    else {
      // fMaxNumberOfSplits < wi_ws
      nw.fN = fMaxNumberOfSplits;
      nw.fW = init_w/fMaxNumberOfSplits;
    }


  } else if (init_w < lowerWeightBound) {
    // Russian roulette
    G4double wi_ws = init_w/survivalWeight;
    G4double p = std::max(wi_ws,1./fMaxNumberOfSplits);
    G4double r = G4UniformRand();
    if (r<p){
      nw.fW = init_w/p;
      nw.fN = 1;
    } else {
      nw.fW = 0;
      nw.fN = 0;
    }
  } 
  // else do nothing
  
  
  return nw;
}

