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
// $Id: G4C11GEMProbability.hh,v 1.1 2002/06/06 17:40:28 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4C11GEMProbability_h
#define G4C11GEMProbability_h 1


#include "G4GEMProbability.hh"


class G4C11GEMProbability : public G4GEMProbability
{
public:
  // Only available constructor
  G4C11GEMProbability();
    
  ~G4C11GEMProbability() {}
private:  
  // Copy constructor
  G4C11GEMProbability(const G4C11GEMProbability &right);
    
  const G4C11GEMProbability & operator=(const G4C11GEMProbability &right);
  G4bool operator==(const G4C11GEMProbability &right) const;
  G4bool operator!=(const G4C11GEMProbability &right) const;
  

private:
  
  virtual G4double CalcBetaParam(const G4Fragment & fragment) const 
  {
    return -GetCoulombBarrier(fragment);
  }
  
  
  // Excitation energy levels 
  G4std::vector<G4double> ExcitEnergies;
  // Spin of excitation energy levels 
  G4std::vector<G4double> ExcitSpins;
  
  G4std::vector<G4double> ExcitLifetimes;
  
};


#endif
