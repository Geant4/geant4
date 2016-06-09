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
// $Id: G4WeightWindowAlgorithm.cc,v 1.8 2006/06/29 18:18:01 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

