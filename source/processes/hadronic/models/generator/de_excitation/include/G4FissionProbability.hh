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
// $Id: G4FissionProbability.hh,v 1.7 2002/12/12 19:17:06 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4FissionProbability_h
#define G4FissionProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4FissionLevelDensityParameter.hh"

class G4FissionProbability : public G4VEmissionProbability
{
public:
  // Default constructor
  G4FissionProbability() {};

  ~G4FissionProbability() {};  

private:  
  // Copy constructor
  G4FissionProbability(const G4FissionProbability &right);

  const G4FissionProbability & operator=(const G4FissionProbability &right);
  G4bool operator==(const G4FissionProbability &right) const;
  G4bool operator!=(const G4FissionProbability &right) const;
  
public:
	G4double EmissionProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy);

private:
	G4EvaporationLevelDensityParameter theEvapLDP;
	G4FissionLevelDensityParameter theFissLDP;


};


#endif
