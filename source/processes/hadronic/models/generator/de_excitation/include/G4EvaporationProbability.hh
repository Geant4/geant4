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
// $Id: G4EvaporationProbability.hh,v 1.4.2.1 2001/06/28 19:13:02 gunter Exp $
// GEANT4 tag $Name:  $
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

  void SetExcitationEnergiesPtr(G4std::vector<G4double> * anExcitationEnergiesPtr) 
  {ExcitationEnergies = anExcitationEnergiesPtr;}

  void SetExcitationSpinsPtr(G4std::vector<G4int> * anExcitationSpinsPtr)
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
  virtual G4double CCoeficient(const G4double aZ) const {return 0.0;};

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const {return 1.0;}
  virtual G4double CalcBetaParam(const G4Fragment & fragment) const {return 1.0;}

  // Data Members

  G4VLevelDensityParameter * theEvapLDPptr;
	
  G4int theA;
  G4int theZ;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;

  // Discrete Excitation Energies 
  G4std::vector<G4double> * ExcitationEnergies;

  //
  G4std::vector<G4int> * ExcitationSpins;

};


#endif
