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
#ifndef G4INUCL_SPECIAL_FUNC_HH
#define G4INUCL_SPECIAL_FUNC_HH

#include "globals.hh"
#include <math.h>
#include "g4std/algorithm"
#include "g4std/vector"


namespace G4InuclSpecialFunctions {

  G4double bindingEnergyExact(G4double A, 
			      G4double Z);

  G4double bindingEnergyKummel(G4double A, 
			       G4double Z);

  G4double bindingEnergy(G4double A, 
			 G4double Z);

  G4double bindingEnergyAsymptotic(G4double A, 
				   G4double Z);

  G4double FermiEnergy(G4double A, 
		       G4double Z, 
		       G4int ntype);
  
  G4std::pair<G4std::vector<G4double>, G4std::vector<G4double> > paraMaker(G4double Z);

  G4std::pair<G4double, G4double> paraMakerTruncated(G4double Z); 

  G4double getAL(G4double A);
 
  G4double csNN(G4double e);

  G4double csPN(G4double e);

  G4double inuclRndm();

  G4double randomGauss(G4double sigma);

  G4double randomPHI();

  G4std::pair<G4double, G4double> randomCOS_SIN();

  G4double nucleiLevelDensity(G4double a);

  G4std::vector<G4double> generateWithFixedTheta(G4double ct, 
					  G4double p);
}
#endif
