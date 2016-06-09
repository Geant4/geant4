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
// $Id$
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
#include "G4VCoulombBarrier.hh"
#include "G4CoulombBarrier.hh"


class G4EvaporationProbability : public G4VEmissionProbability
{
public:
  // Only available constructor
  G4EvaporationProbability(G4int anA, G4int aZ, G4double aGamma,
			   G4VCoulombBarrier * aCoulombBarrier); 

  virtual ~G4EvaporationProbability();

  inline G4int GetZ(void) const { return theZ; }
	
  inline G4int GetA(void) const { return theA;} 

protected:
  
  // Default constructor
  G4EvaporationProbability();

private:
  // Copy constructor
  G4EvaporationProbability(const G4EvaporationProbability &right);

  const G4EvaporationProbability & operator=(const G4EvaporationProbability &right);
  G4bool operator==(const G4EvaporationProbability &right) const;
  G4bool operator!=(const G4EvaporationProbability &right) const;
  
public:

  G4double ProbabilityDistributionFunction( const G4Fragment & aFragment, G4double K);

  G4double EmissionProbability(const G4Fragment & fragment, G4double anEnergy);

private:

  G4double CalculateProbability(const G4Fragment & fragment, G4double MaximalKineticEnergy );

  G4double IntegrateEmissionProbability(const G4Fragment & aFragment, 
					const G4double & Low, const G4double & Up );

protected:

 virtual G4double CrossSection( const  G4Fragment & fragment, G4double K )= 0;  

 virtual G4double CalcAlphaParam(const G4Fragment & fragment)=0 ;
 
 virtual G4double CalcBetaParam(const G4Fragment & fragment)=0 ;

private:

  // Data Members	
  G4int theA;
  G4int theZ;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;

  //The Coulomb Barrier
  G4VCoulombBarrier * theCoulombBarrierptr;

};

#endif
