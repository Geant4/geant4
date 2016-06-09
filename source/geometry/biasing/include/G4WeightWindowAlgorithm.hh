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
// $Id: G4WeightWindowAlgorithm.hh,v 1.6 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowAlgorithm
// 
// Class description:
// Implementation of a weight window algorithm. 
// The arguments in the constructor configure the algorithm:
//   - upperLimitFaktor: the factor defining the upper weight limit 
//      W_u = upperLimitFaktor * W_l (W_l lower weight bound)
//   - survivalFaktor: used in calculating  the survival weight
//      W_s = survivalFaktor * W_l
//   - maxNumberOfSplits: the maximal number of splits allowed
//     to be created in one go, and the reciprocal of the minimal
//     survival probability in case of Russian roulette
//
// In case of upperLimitFaktor=survivalFaktor=1 the algorithm
// becomes the expected weight algorithm of the importance sampling 
// technique.
//
// see also G4VWeightWindowAlgorithm.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightWindowAlgorithm_hh
#define G4WeightWindowAlgorithm_hh G4WeightWindowAlgorithm_hh 

#include "G4Nsplit_Weight.hh"
#include "G4VWeightWindowAlgorithm.hh"

class G4WeightWindowAlgorithm : public G4VWeightWindowAlgorithm
{

public:  // with description
  
  G4WeightWindowAlgorithm(G4double upperLimitFaktor = 5,
			  G4double survivalFaktor = 3,
			  G4int maxNumberOfSplits = 5);
  
  virtual ~G4WeightWindowAlgorithm();

  virtual G4Nsplit_Weight Calculate(G4double init_w,
				    G4double lowerWeightBound) const;
    // calculate number of tracks and their weight according
    // to the initial track weight and the lower energy bound

private:
  G4double fUpperLimitFaktor;
  G4double fSurvivalFaktor;
  G4int fMaxNumberOfSplits;
};

#endif
