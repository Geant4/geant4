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
// $Id: G4VEmissionProbability.hh,v 1.2 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4VEmissionProbability_h
#define G4VEmissionProbability_h 1


#include "globals.hh"
#include "G4Fragment.hh"

class G4VEmissionProbability 
{
public:
  G4VEmissionProbability() {};
  virtual ~G4VEmissionProbability() {};

private:  
  G4VEmissionProbability(const G4VEmissionProbability &right);

  const G4VEmissionProbability & operator=(const G4VEmissionProbability &right);
  G4bool operator==(const G4VEmissionProbability &right) const;
  G4bool operator!=(const G4VEmissionProbability &right) const;
  
public:
  virtual G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy) = 0;

};


#endif
