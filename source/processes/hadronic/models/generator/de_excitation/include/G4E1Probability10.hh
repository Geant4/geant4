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
//

#ifndef G4E1Probability10_hh
#define G4E1Probability10_hh

#include "globals.hh"
#include "G4VEmissionProbability.hh"
#include "G4Fragment.hh"
#include "G4VLevelDensityParameter.hh"

class G4E1Probability10 : public G4VEmissionProbability
{

public:

  G4E1Probability10() {};

  ~G4E1Probability10();

  G4double EmissionProbability(const G4Fragment& frag, const G4double excite);
  G4double EmissionProbDensity(const G4Fragment& frag, const G4double ePhoton);

private:

  // G4E1Probability10() {};

  G4E1Probability10(const G4E1Probability10& right);

  const G4E1Probability10& operator=(const G4E1Probability10& right);
  G4bool operator==(const G4E1Probability10& right) const;
  G4bool operator!=(const G4E1Probability10& right) const;

  // Integrator (simple Gaussian quadrature)

  G4double EmissionIntegration(const G4Fragment& frag, const G4double excite,
                               const G4double lowLim, const G4double upLim,
                               const G4int numIters);

  // G4VLevelDensityParameter* _levelDensity; // Don't need this

};

#endif
