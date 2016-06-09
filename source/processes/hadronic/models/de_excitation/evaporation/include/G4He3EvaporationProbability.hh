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
// $Id: G4He3EvaporationProbability.hh,v 1.1 2003/08/26 18:30:46 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4He3EvaporationProbability_h
#define G4He3EvaporationProbability_h 1


#include "G4EvaporationProbability.hh"


class G4He3EvaporationProbability : public G4EvaporationProbability
{
public:
  // Only available constructor
  G4He3EvaporationProbability();

  ~G4He3EvaporationProbability() {}
private:  
  // Copy constructor
  G4He3EvaporationProbability(const G4He3EvaporationProbability &right);

  const G4He3EvaporationProbability & operator=(const G4He3EvaporationProbability &right);
  G4bool operator==(const G4He3EvaporationProbability &right) const;
  G4bool operator!=(const G4He3EvaporationProbability &right) const;
  

private:

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
  { return 1.0 + CCoeficient(static_cast<G4double>(fragment.GetZ()-GetZ()));}
	
  virtual G4double CalcBetaParam(const G4Fragment & ) const 
  { return 0.0; }

	
  G4double CCoeficient(const G4double aZ) const;

  // Excitation energy levels 
  std::vector<G4double> ExcitEnergies;
  // Spin of excitation energy levels 
  std::vector<G4int> ExcitSpins;
};


#endif
