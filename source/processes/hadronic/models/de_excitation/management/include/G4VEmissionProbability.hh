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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICB is set true, by default is false) 
//
// V.Ivanchenko general clean-up since 2010
//

#ifndef G4VEmissionProbability_h
#define G4VEmissionProbability_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4NuclearLevelData;
class G4Pow;

class G4VEmissionProbability 
{
public:

  explicit G4VEmissionProbability(G4int Z, G4int A);

  virtual ~G4VEmissionProbability() = default;

  void Initialise();

  virtual G4double EmissionProbability(const G4Fragment & fragment, 
  				       G4double anEnergy);

  virtual G4double ComputeProbability(G4double anEnergy, G4double CB);

  inline G4int GetZ(void) const { return theZ; }
	
  inline G4int GetA(void) const { return theA; }

  // Z, A, rmass are residual parameters
  // fmass is SCM mass of decaying nucleus
  // exc is an excitation of emitted fragment
  inline void SetDecayKinematics(G4int rZ, G4int rA, G4double rmass, 
                                 G4double fmass);

  inline G4double GetRecoilExcitation() const { return fExcRes; };

  inline void SetEvapExcitation(G4double exc) { fExc = exc; };

  inline G4double GetProbability() const { return pProbability; };

  inline void ResetProbability() { pProbability = 0.0; };

  // this method may be called only if the probability is computed
  // for given initial fragment and decay channel
  G4double SampleEnergy();

  G4VEmissionProbability(const G4VEmissionProbability &right) = delete;
  const G4VEmissionProbability & operator=
  (const G4VEmissionProbability &right) = delete;
  G4bool operator==(const G4VEmissionProbability &right) const = delete;
  G4bool operator!=(const G4VEmissionProbability &right) const = delete;

protected:

  void ResetIntegrator(size_t nbin, G4double de, G4double eps);

  G4double IntegrateProbability(G4double elow, G4double ehigh, G4double CB);

  G4NuclearLevelData* pNuclearLevelData;
  G4Pow* pG4pow;

  G4int OPTxs;
  G4int pVerbose;
  G4int theZ;
  G4int theA;
  G4int resZ = 0;
  G4int resA = 0;

  G4double pMass = 0.0; // initial fragment
  G4double pEvapMass = 0.0;
  G4double pResMass = 0.0;
  G4double pProbability = 0.0;
  G4double pTolerance = 0.0;

private:

  G4double FindRecoilExcitation(const G4double e);

  G4double fExc = 0.0;
  G4double fExcRes = 0.0;

  G4double fE1 = 0.0;
  G4double fE2 = 0.0; 
  G4double fP2 = 0.0; 

  G4double emin = 0.0;
  G4double emax = 0.0;
  G4double eCoulomb = 0.0;
  G4double accuracy = 0.005;
  G4double probmax = 0.0;
  G4double elimit;

  G4bool fFD = false;
};

inline void 
G4VEmissionProbability::SetDecayKinematics(G4int rZ, G4int rA, G4double rmass, 
                                           G4double fmass)
{
  resZ = rZ;
  resA = rA;
  pMass = fmass;
  pResMass = rmass;
}

#endif
