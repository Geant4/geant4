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
//
// $Id: G4WeightWindowAlgorithm.cc 66356 2012-12-18 09:02:32Z gcosmo $
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

    //TB    G4double wi_ws = init_w/survivalWeight;
    //TB
    G4double temp_wi_ws = init_w/upperWeight;
    G4int split_i = static_cast<int>(temp_wi_ws);
    if(split_i != temp_wi_ws) split_i++;
    G4double wi_ws = init_w/split_i;

    //TB
//TB    G4int int_wi_ws = static_cast<int>(wi_ws);

    // values in case integer mode or in csae of double
    // mode and the lower number of splits will be diced
    //TB    nw.fN = int_wi_ws;
    //TB    nw.fW = survivalWeight;	
    nw.fN = split_i;
    nw.fW = wi_ws;	

//TB     if (wi_ws <= fMaxNumberOfSplits) {
//TB       if (wi_ws > int_wi_ws) {
//TB 	// double mode
//TB 	G4double p2 =  wi_ws - int_wi_ws;
//TB 	G4double r = G4UniformRand();
//TB 	if (r<p2) {
//TB 	  nw.fN = int_wi_ws + 1;
//TB 	}
//TB       }
//TB     }
//TB     else {
//TB       // fMaxNumberOfSplits < wi_ws
//TB       nw.fN = fMaxNumberOfSplits;
//TB       nw.fW = init_w/fMaxNumberOfSplits;
//TB     }


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

