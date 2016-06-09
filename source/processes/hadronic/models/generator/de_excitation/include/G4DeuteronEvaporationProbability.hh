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
// $Id: G4DeuteronEvaporationProbability.hh,v 1.5 2002/12/12 19:17:03 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4DeuteronEvaporationProbability_h
#define G4DeuteronEvaporationProbability_h 1


#include "G4EvaporationProbability.hh"


class G4DeuteronEvaporationProbability : public G4EvaporationProbability
{
public:
  // Only available constructor
  G4DeuteronEvaporationProbability();

  ~G4DeuteronEvaporationProbability() {}
private:  
  // Copy constructor
  G4DeuteronEvaporationProbability(const G4DeuteronEvaporationProbability &right);

  const G4DeuteronEvaporationProbability & operator=(const G4DeuteronEvaporationProbability &right);
  G4bool operator==(const G4DeuteronEvaporationProbability &right) const;
  G4bool operator!=(const G4DeuteronEvaporationProbability &right) const;
  

private:

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
  { return 1.0 + CCoeficient(G4double(fragment.GetZ()-GetZ()));}
	
  virtual G4double CalcBetaParam(const G4Fragment & fragment) const 
  { return 0.0; }

	
  G4double CCoeficient(const G4double aZ) const;

  // Excitation energy levels 
  G4std::vector<G4double> ExcitEnergies;
  // Spin of excitation energy levels 
  G4std::vector<G4int> ExcitSpins;

};
#endif
