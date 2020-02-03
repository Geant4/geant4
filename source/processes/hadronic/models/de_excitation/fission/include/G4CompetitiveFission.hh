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
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)

#ifndef G4CompetitiveFission_h
#define G4CompetitiveFission_h 1

#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4VEmissionProbability.hh"
#include "G4FissionParameters.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4Exp.hh"

class G4VFissionBarrier;
class G4VEmissionProbability;
class G4VLevelDensityParameter;
class G4PairingCorrection;

class G4CompetitiveFission : public G4VEvaporationChannel
{
public:
  
  explicit G4CompetitiveFission();
  ~G4CompetitiveFission() override;

  G4Fragment* EmittedFragment(G4Fragment* theNucleus) override;

  G4double GetEmissionProbability(G4Fragment* theNucleus) override;

  void SetFissionBarrier(G4VFissionBarrier * aBarrier);

  void SetEmissionStrategy(G4VEmissionProbability * aFissionProb);

  void SetLevelDensityParameter(G4VLevelDensityParameter * aLevelDensity);

  inline G4double GetFissionBarrier(void) const;

  inline G4double GetLevelDensityParameter(void) const;

  inline G4double GetMaximalKineticEnergy(void) const;

private:

  // Sample AtomicNumber of Fission products
  G4int FissionAtomicNumber(G4int A);

  G4double MassDistribution(G4double x, G4int A);

  // Sample Charge of fission products
  G4int FissionCharge(G4int A, G4int Z, G4double Af);

  // Sample Kinetic energy of fission products
  G4double FissionKineticEnergy(G4int A, G4int Z,
				G4int Af1, G4int Zf1,
				G4int Af2, G4int Zf2,
				G4double U, G4double Tmax);
    
  inline G4double Ratio(G4double A, G4double A11, 
                        G4double B1, G4double A00) const;

  inline G4double SymmetricRatio(G4int A, G4double A11) const;

  inline G4double AsymmetricRatio(G4int A, G4double A11) const;

  inline G4double LocalExp(G4double x) const;

  G4CompetitiveFission(const G4CompetitiveFission &right);
  const G4CompetitiveFission & operator=(const G4CompetitiveFission &right);
  G4bool operator==(const G4CompetitiveFission &right) const;
  G4bool operator!=(const G4CompetitiveFission &right) const;

  // Maximal Kinetic Energy that can be carried by fragment
  G4double maxKineticEnergy;
  G4double fissionBarrier;
  G4double fissionProbability;

  // For Fission barrier
  G4VFissionBarrier* theFissionBarrierPtr;

  // For Fission probability emission
  G4VEmissionProbability* theFissionProbabilityPtr;

  // For Level Density calculation
  G4VLevelDensityParameter* theLevelDensityPtr;
  G4PairingCorrection* pairingCorrection;

  G4bool myOwnFissionProbability;
  G4bool myOwnFissionBarrier;
  G4bool myOwnLevelDensity;

  G4FissionParameters theParam;

};

inline G4double G4CompetitiveFission::GetFissionBarrier(void) const 
{ 
  return fissionBarrier; 
}

inline G4double G4CompetitiveFission::GetMaximalKineticEnergy(void) const 
{ 
  return maxKineticEnergy; 
}

inline
G4double G4CompetitiveFission::Ratio(G4double A, G4double A11,
				     G4double B1, G4double A00) const
{
  G4double res;
  if (A11 >= A*0.5 && A11 <= (A00+10.0)) {
    G4double x = (A11-A00)/A;
    res = 1.0 - B1*x*x;
  } else {
    G4double x = 10.0/A;
    res = 1.0 - B1*x*x - 2.0*x*B1*(A11-A00-10.0)/A;
  }
  return res;
}

inline
G4double G4CompetitiveFission::AsymmetricRatio(G4int A, G4double A11) const
{
  return Ratio(G4double(A),A11,23.5,134.0);
}

inline
G4double G4CompetitiveFission::SymmetricRatio(G4int A, G4double A11) const
{
  G4double A0 = G4double(A);
  return Ratio(A0,A11,5.32,A0*0.5);
}

inline G4double G4CompetitiveFission::LocalExp(G4double x) const
{
  return (std::abs(x) < 8.) ? G4Exp(-0.5*x*x) : 0.0;
}

#endif


