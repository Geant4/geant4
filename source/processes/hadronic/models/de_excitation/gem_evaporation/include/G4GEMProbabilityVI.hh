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
// GEM de-excitation model
// by V. Ivanchenko (July 2019)
//
#ifndef G4GEMProbabilityVI_h
#define G4GEMProbabilityVI_h 1

#include "G4VEmissionProbability.hh"

class G4LevelManager;

class G4GEMProbabilityVI : public G4VEmissionProbability
{
public:

  explicit G4GEMProbabilityVI(G4int anA, G4int aZ, const G4LevelManager*); 

  ~G4GEMProbabilityVI() override = default;

  G4double TotalProbability(const G4Fragment&,
                            const G4double tmin, const G4double tmax, 
                            const G4double CB, const G4double exEnergy,
                            const G4double exEvap);

  // compute probability for evaporated fragment in ground state
  G4double ComputeProbability(G4double ekin, G4double CB) override;

  G4double SampleEnergy(const G4double tmin, const G4double tmax, 
			const G4double CB, const G4double exEnergy,
			const G4double exEvap);

  G4GEMProbabilityVI(const G4GEMProbabilityVI& right) = delete;
  const G4GEMProbabilityVI & operator=(const G4GEMProbabilityVI& right) = delete;
  G4bool operator==(const G4GEMProbabilityVI& right) const = delete;
  G4bool operator!=(const G4GEMProbabilityVI& right) const = delete;

private:

  const G4LevelManager* lManager;

  G4int fragA;
  G4int fragZ;

  G4double bCoulomb;
  G4double resA13;
  G4double U, delta0, delta1, a0, a1;
  G4double alphaP, betaP;
  G4double Umax, A13;
  //  G4double levelDensity, levelDensity1;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;
  G4double coeff;
  G4double pcoeff;

  G4double probmax;

  G4bool isExcited;

  //static const G4double ws[NPOINTSGEM];
  //static const G4double xs[NPOINTSGEM];

};

#endif
