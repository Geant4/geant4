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
#ifndef G4INUCL_SPECIAL_FUNC_HH
#define G4INUCL_SPECIAL_FUNC_HH

#include "globals.hh"
#include <cmath>
#include <algorithm>
#include <vector>
#include "G4CascadeMomentum.hh"

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
  
  std::pair<std::vector<G4double>, std::vector<G4double> > paraMaker(G4double Z);

  std::pair<G4double, G4double> paraMakerTruncated(G4double Z); 

  G4double getAL(G4double A);
 
  G4double csNN(G4double e);

  G4double csPN(G4double e);

  G4double inuclRndm();

  G4double randomGauss(G4double sigma);

  G4double randomPHI();

  std::pair<G4double, G4double> randomCOS_SIN();

  G4double nucleiLevelDensity(G4double a);

  G4CascadeMomentum generateWithFixedTheta(G4double ct, 
					  G4double p);
}
#endif
