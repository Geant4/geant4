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
// $Id: G4WeightWindowAlgorithm.hh,v 1.5 2003-08-15 15:34:45 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowAlgorithm
// 
// Class description:
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

private:
  G4double fUpperLimitFaktor;
  G4double fSurvivalFaktor;
  G4int fMaxNumberOfSplits;
};

#endif
