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
// $Id: G4EvaporationProbability.hh,v 1.3 2006/06/29 20:09:57 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4EvaporationProbability_h
#define G4EvaporationProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"


class G4EvaporationProbability : public G4VEmissionProbability
{
public:
  // Only available constructor
  G4EvaporationProbability(const G4int anA, const G4int aZ, const G4double aGamma) : 
    theA(anA),
    theZ(aZ),
    Gamma(aGamma) 
  {
    theEvapLDPptr = new G4EvaporationLevelDensityParameter;
  }

  ~G4EvaporationProbability() 
  {
    if (theEvapLDPptr != 0) delete theEvapLDPptr;
  }


	
  G4double GetZ(void) const { return theZ; }
	
  G4double GetA(void) const { return theA;} 

protected:

  void SetExcitationEnergiesPtr(std::vector<G4double> * anExcitationEnergiesPtr) 
  {ExcitationEnergies = anExcitationEnergiesPtr;}

  void SetExcitationSpinsPtr(std::vector<G4int> * anExcitationSpinsPtr)
  {ExcitationSpins = anExcitationSpinsPtr;}

  
  // Default constructor
  G4EvaporationProbability() {}
private:
  // Copy constructor
  G4EvaporationProbability(const G4EvaporationProbability &right);

  const G4EvaporationProbability & operator=(const G4EvaporationProbability &right);
  G4bool operator==(const G4EvaporationProbability &right) const;
  G4bool operator!=(const G4EvaporationProbability &right) const;
  
public:
  G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy);

private:

  G4double CalcProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy);
  virtual G4double CCoeficient(const G4double ) const {return 0.0;};

  virtual G4double CalcAlphaParam(const G4Fragment & ) const {return 1.0;}
  virtual G4double CalcBetaParam(const G4Fragment & ) const {return 1.0;}

  // Data Members

  G4VLevelDensityParameter * theEvapLDPptr;
	
  G4int theA;
  G4int theZ;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;

  // Discrete Excitation Energies 
  std::vector<G4double> * ExcitationEnergies;

  //
  std::vector<G4int> * ExcitationSpins;

};


#endif
