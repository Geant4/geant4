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
#include <vector>

class G4NuclearLevelData;
class G4Pow;

class G4VEmissionProbability 
{
public:

  explicit G4VEmissionProbability(G4int Z, G4int A);
  virtual ~G4VEmissionProbability();

  void Initialise();

  virtual G4double EmissionProbability(const G4Fragment & fragment, 
  				       G4double anEnergy);

  virtual G4double ComputeProbability(G4double anEnergy, G4double CB);

  inline G4int GetZ(void) const { return theZ; }
	
  inline G4int GetA(void) const { return theA; }

  // Z, A, rmass are residual parameters
  // fmass is SCM mass of decaying nucleus
  // exc is an excitation of emitted fragment
  inline void SetDecayKinematics(G4int Z, G4int A, G4double rmass, 
                                 G4double fmass);

  inline G4double GetRecoilExcitation() const { return fExcRes; };

  inline void SetEvapExcitation(G4double exc) { fExc = exc; };

  inline G4double GetProbability() const { return pProbability; };

  inline void ResetProbability() { pProbability = 0.0; };

  // this method may be called only if the probability is computed
  // for given initial fragment and decay channel
  G4double SampleEnergy();

protected:

  void ResetIntegrator(size_t nbin, G4double de, G4double eps);

  G4double IntegrateProbability(G4double elow, G4double ehigh, G4double CB);

  G4int OPTxs;
  G4int pVerbose;
  G4int theZ;
  G4int theA;
  G4int resZ;
  G4int resA;

  G4double pMass; // initial fragment
  G4double pEvapMass;
  G4double pResMass;
  G4double pProbability;

  G4NuclearLevelData* pNuclearLevelData;
  G4Pow* pG4pow;

private:

  G4double FindRecoilExcitation(G4double e);

  G4VEmissionProbability(const G4VEmissionProbability &right);
  const G4VEmissionProbability & operator=
  (const G4VEmissionProbability &right);
  G4bool operator==(const G4VEmissionProbability &right) const;
  G4bool operator!=(const G4VEmissionProbability &right) const;

  size_t length;
  size_t nbin;

  G4double fExc;
  G4double fExcRes;

  G4double emin;
  G4double emax;
  G4double elimit;
  G4double eCoulomb;
  G4double accuracy;
  G4double probmax;

  G4bool   fFD;
};

inline void 
G4VEmissionProbability::SetDecayKinematics(G4int Z, G4int A, G4double rmass, 
                                           G4double fmass)
{
  resZ = Z;
  resA = A;
  pMass = fmass;
  pResMass = rmass;
}

#endif
