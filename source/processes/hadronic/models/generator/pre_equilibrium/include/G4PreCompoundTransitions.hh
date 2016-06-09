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
// $Id: G4PreCompoundTransitions.hh,v 1.10 2003/05/30 14:05:57 hpw Exp $
// GEANT4 tag $Name: geant4-05-02 $
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

#include "G4VPreCompoundTransitions.hh"

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

class G4PreCompoundTransitions : public G4VPreCompoundTransitions
{
public:

  // Calculates transition probabilities with Delta N = +2 (Trans1) -2 (Trans2) and 0 (Trans3)

  G4PreCompoundTransitions() : TransitionProb1(0.0), TransitionProb2(0.0), TransitionProb3(0.0) {}

  virtual ~G4PreCompoundTransitions() {}

private:
  
  G4PreCompoundTransitions(const G4PreCompoundTransitions &) : G4VPreCompoundTransitions() {};
  
  const G4PreCompoundTransitions& operator=(const G4PreCompoundTransitions &right);

  G4bool operator==(const G4PreCompoundTransitions &right) const;
  
  G4bool operator!=(const G4PreCompoundTransitions &right) const;

public:
  
  virtual G4double CalculateProbability(const G4Fragment & aFragment);
  
  virtual G4Fragment PerformTransition(const G4Fragment & aFragment);
  
private:
	
  G4double TransitionProb1;
  G4double TransitionProb2;
  G4double TransitionProb3;

};

#endif
