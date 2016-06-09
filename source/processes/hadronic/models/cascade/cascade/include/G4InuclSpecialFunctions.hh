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
// $Id: G4InuclSpecialFunctions.hh,v 1.18 2010-09-14 17:51:36 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100319  M. Kelsey -- Add optional mass argument to generateWithFixedTheta;
//		define new generateWithRandomAngles, encapsulating code; define
//		cbrt() cube-root function (in math.h, but not in <math>!)
// 20100412  M. Kelsey -- Modify paraMaker[Truncated] to take buffer as argument
// 20100914  M. Kelsey -- Migrate to integer A and Z.  Discard unused binding
//		energy functions

#ifndef G4INUCL_SPECIAL_FUNC_HH
#define G4INUCL_SPECIAL_FUNC_HH

#include "globals.hh"
#include <algorithm>
#include <vector>
#include "G4LorentzVector.hh"

namespace G4InuclSpecialFunctions {
  G4double bindingEnergy(G4int A, G4int Z);

  // NOTE:  Used only by G4Fissioner
  G4double bindingEnergyAsymptotic(G4int A, G4int Z);

  G4double FermiEnergy(G4int A, G4int Z, G4int ntype);
  
  // NOTE:  Passing Z as double here, to be used as interpolation argument
  void paraMaker(G4double Z, std::pair<std::vector<G4double>, std::vector<G4double> >& parms);

  void paraMakerTruncated(G4double Z, std::pair<G4double, G4double>& parms); 

  G4double getAL(G4int A);
 
  G4double csNN(G4double e);

  G4double csPN(G4double e);

  G4double G4cbrt(G4double x);	// Can't use "cbrt" name, clashes with <math.h>

  G4double inuclRndm();

  G4double randomGauss(G4double sigma);

  G4double randomPHI();

  std::pair<G4double, G4double> randomCOS_SIN();

  G4double nucleiLevelDensity(G4int A);

  // Optional mass argument will be used to fill G4LorentzVector correctly
  G4LorentzVector generateWithFixedTheta(G4double ct, G4double p,
					 G4double m=0.);

  G4LorentzVector generateWithRandomAngles(G4double p, G4double m=0.);
}


#endif
