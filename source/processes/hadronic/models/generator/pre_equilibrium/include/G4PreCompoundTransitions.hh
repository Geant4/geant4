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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PreCompoundTransitions.hh,v 1.7 2001/08/01 17:08:29 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// by V. Lara

#ifndef G4PreCompoundTransitions_h
#define G4PreCompoundTransitions_h 1

// Compute transition probailities:
// TransitionProb1 => probability of transition with  \Delta N = +2
//                    number of excitons will be increased on 2
// TransitionProb2 => probability of transition with  \Delta N = -2
//                    number of excitons will be decreased on 2
// TransitionProb3 => probability of transition with  \Delta N = 0
//                    number of excitons will be the same

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

class G4PreCompoundTransitions
{
public:

  // Calculates transition probabilities with Delta N = +2 (Trans1) -2 (Trans2) and 0 (Trans3)
  G4PreCompoundTransitions(const G4Fragment & aFragment);

  ~G4PreCompoundTransitions() {}

private:
  G4PreCompoundTransitions() {}
  
  G4PreCompoundTransitions(const G4PreCompoundTransitions &right) {};
  
  const G4PreCompoundTransitions& operator=(const G4PreCompoundTransitions &right);

  G4bool operator==(const G4PreCompoundTransitions &right) const;
  
  G4bool operator!=(const G4PreCompoundTransitions &right) const;

public:
  G4double GetTotalProbability(void)
  { return TransitionProb1+TransitionProb2+TransitionProb3; }
  
  G4Fragment PerformTransition(const G4Fragment & aFragment);
  
private:
	
  G4double TransitionProb1;
  G4double TransitionProb2;
  G4double TransitionProb3;

};

#endif
