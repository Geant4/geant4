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
// $Id: G4EvaporationProbability.hh 103162 2017-03-20 09:40:58Z gcosmo $
//
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
#ifndef G4EvaporationProbability_h
#define G4EvaporationProbability_h 1

#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"

class G4VCoulombBarrier;
class G4NuclearLevelData;

class G4EvaporationProbability : public G4VEmissionProbability
{
public:

  explicit G4EvaporationProbability(G4int anA, G4int aZ, 
                                    G4double aGamma, 
                                    G4VCoulombBarrier *); 

  virtual ~G4EvaporationProbability();

  // not used for evaporation
  virtual G4double EmissionProbability(const G4Fragment& fragment,
				       G4double maxKineticEnergy);

  // general method used for evaporation
  G4double TotalProbability(const G4Fragment& fragment,
			    G4double minKineticEnergy,
			    G4double maxKineticEnergy,
			    G4double CoulombBarrier = 0.0);

  // main method to compute full probability for OPTx > 2
  virtual G4double ComputeProbability(G4double K, G4double kBarrier);

  // Samples fragment kinetic energy and excitation energy 
  // of the residual nucleaus
  G4double SampleKineticEnergy(G4double minKineticEnergy,
			       G4double maxKineticEnergy,
			       G4double CoulombBarrier = 0.0);

protected:

  virtual G4double CalcAlphaParam(const G4Fragment & fragment)=0 ;
 
  virtual G4double CalcBetaParam(const G4Fragment & fragment)=0 ;

private:

  G4double CrossSection(G4double K, G4double CoulombBarrier);  

  // Copy constructor
  G4EvaporationProbability(const G4EvaporationProbability &right) = delete;

  const G4EvaporationProbability & operator=
  (const G4EvaporationProbability &right) = delete;
  G4bool operator==(const G4EvaporationProbability &right) const = delete;
  G4bool operator!=(const G4EvaporationProbability &right) const = delete;

  G4NuclearLevelData* fLevelData;

  G4int fragA;
  G4int fragZ;
  G4int resA;
  G4int resZ;
  G4int index;

  G4double resA13;
  G4double muu;
  G4double partMass;
  G4double resMass;
  G4double Mass;
  G4double U, delta0, delta1, a0;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;
  G4double pcoeff;

  G4double probmax;
};

#endif
