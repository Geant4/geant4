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
// $Id: G4WeightWindowAlgorithm.cc,v 1.4 2002-10-14 12:36:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WeightWindowAlgorithm.cc
//
// ----------------------------------------------------------------------
#include "G4WeightWindowAlgorithm.hh"
#include "Randomize.hh"


G4WeightWindowAlgorithm::G4WeightWindowAlgorithm() :
  fUpper(1),
  fLower(1)
{}

G4WeightWindowAlgorithm::~G4WeightWindowAlgorithm()
{}


void G4WeightWindowAlgorithm::SetUpperLimit(G4double Upper){
  fUpper = Upper;
}
  
void G4WeightWindowAlgorithm::SetLowerLimit(G4double Lower){
  fLower = Lower;
}


G4Nsplit_Weight 
G4WeightWindowAlgorithm::Calculate(G4double init_w, 
				   G4double importance) const {


  G4Nsplit_Weight nw(1, init_w);
  G4double iw =  importance * init_w;
  if (iw>fUpper) {
    // f is the factor by which the weight is greater 
    // than allowed by fUpper
    // it is almost the number of coppies to be produced
    G4double f = iw / fUpper;
    // calculate new weight
    nw.fW/=f; 

    // calculate the number of coppies
    nw.fN = static_cast<G4int>(f);
    // correct the number of coppies in case f is not an integer
    if (static_cast<G4double>(nw.fN) != f) {
      // probabillity p for splitting into nw.fN+1 particles
      G4double p = f - nw.fN;
      // get a random number out of [0,1)
      G4double r = G4UniformRand();
      if (r<p) {
	nw.fN++;
      } 
    }

  }
  else if (iw < fLower) {
    // play russian roulett
    G4double p = iw / fLower; // survival prob.
    G4double r = G4UniformRand();
    if (r>=p) {
      // kill track
      nw.fN = 0;
      nw.fW = 0;
    } 
    else {
      nw.fW*=1/p; // to be consistant with the survival prob.
    }      
  }
  return nw;
}

