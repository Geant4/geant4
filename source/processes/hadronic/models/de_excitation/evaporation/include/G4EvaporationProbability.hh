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
// $Id: G4EvaporationProbability.hh 90273 2015-05-22 10:20:32Z gcosmo $
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

class G4EvaporationProbability : public G4VEmissionProbability
{
public:

  G4EvaporationProbability(G4int anA, G4int aZ, G4double aGamma,
			   G4VCoulombBarrier *); 

  virtual ~G4EvaporationProbability();

  inline G4int GetZ(void) const { return theZ; }
	
  inline G4int GetA(void) const { return theA; }

  // obsolete method
  G4double EmissionProbability(const G4Fragment& fragment,
			       G4double maxKineticEnergy);

  G4double TotalProbability(const G4Fragment& fragment,
			    G4double minKineticEnergy,
			    G4double maxKineticEnergy);

  G4double ProbabilityDistributionFunction(G4double K);

  // Samples fragment kinetic energy.
  G4double SampleKineticEnergy(G4double minKineticEnergy,
			       G4double maxKineticEnergy);

protected:

  virtual G4double CalcAlphaParam(const G4Fragment & fragment)=0 ;
 
  virtual G4double CalcBetaParam(const G4Fragment & fragment)=0 ;

private:

  G4double IntegrateEmissionProbability(G4double low, G4double up);

  G4double CrossSection(G4double K);  

  // Copy constructor
  G4EvaporationProbability(const G4EvaporationProbability &right);

  const G4EvaporationProbability & operator=(const G4EvaporationProbability &right);
  G4bool operator==(const G4EvaporationProbability &right) const;
  G4bool operator!=(const G4EvaporationProbability &right) const;

  G4int theA;
  G4int theZ;
  G4int fragA;
  G4int fragZ;
  G4int resA;
  G4int resZ;
  G4int index;
  G4int nbins;

  G4double resA13;
  G4double muu;
  G4double partMass;
  G4double resMass;
  G4double fragMass;
  G4double U, delta0, delta1, a0;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;

  G4double probability[11];
};

#endif
