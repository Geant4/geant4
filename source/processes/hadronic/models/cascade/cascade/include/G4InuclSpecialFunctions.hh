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
// G4InuclSpecialFunctions namespace
//
// Original Author: Aatos Heikkinen, 2002
// --------------------------------------------------------------------
#ifndef G4INUCL_SPECIAL_FUNC_HH
#define G4INUCL_SPECIAL_FUNC_HH

#include "globals.hh"
#include "G4LorentzVector.hh"
#include <algorithm>
#include <vector>

template <int N> class G4CascadeInterpolator;

namespace G4InuclSpecialFunctions
{
  G4double bindingEnergy(G4int A, G4int Z);

  // NOTE:  Used only by G4Fissioner
  G4double bindingEnergyAsymptotic(G4int A, G4int Z);

  G4double FermiEnergy(G4int A, G4int Z, G4int ntype);
  
  G4double getAL(G4int A);
 
  G4double csNN(G4double e);

  G4double csPN(G4double e);

  G4double G4cbrt(G4double x);	// Use std::cbrt from <cmath>
  G4double G4cbrt(G4int n);	// Use G4Pow::powN() here for speedup

  G4double randomInuclPowers(G4double ekin, 	// Power series in Ekin, S
			     const G4double (&coeff)[4][4]);

  G4double randomGauss(G4double sigma);		// Gaussian distribution

  G4double randomPHI();

  std::pair<G4double, G4double> randomCOS_SIN();

  G4double nucleiLevelDensity(G4int A);

  // Optional mass argument will be used to fill G4LorentzVector correctly
  G4LorentzVector generateWithFixedTheta(G4double ct, G4double p,
					 G4double mass=0.);

  G4LorentzVector generateWithRandomAngles(G4double p, G4double mass=0.);

  // This is an object to be instantiated by client code
  class paraMaker
  {
    public:

      paraMaker(G4int verbose=0);
      ~paraMaker();

      paraMaker(const paraMaker& right) = delete;
      paraMaker& operator=(const paraMaker& right) = delete;

      // NOTE:  Passing Z as double here, to be used as interpolation argument
      void getParams(G4double Z, std::pair<std::vector<G4double>,
                                 std::vector<G4double> >& parms);

      void getTruncated(G4double Z, std::pair<G4double, G4double>& parms); 

    private:

      G4int verboseLevel;
      G4CascadeInterpolator<5>* interp;
  };
}


#endif
