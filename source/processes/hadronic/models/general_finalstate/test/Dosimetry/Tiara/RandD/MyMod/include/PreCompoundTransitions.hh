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
// $Id: PreCompoundTransitions.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara

#ifndef PreCompoundTransitions_h
#define PreCompoundTransitions_h 1

// Compute transition probailities:
// TransitionProb1 => probability of transition with  \Delta N = +2
//                    number of excitons will be increased on 2
// TransitionProb2 => probability of transition with  \Delta N = -2
//                    number of excitons will be decreased on 2
// TransitionProb3 => probability of transition with  \Delta N = 0
//                    number of excitons will be the same

#include "globals.hh"
#include "G4Fragment.hh"
#include "PreCompoundParameters.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

class PreCompoundTransitions
{
public:

  // Calculates transition probabilities with Delta N = +2 (Trans1) -2 (Trans2) and 0 (Trans3)
  PreCompoundTransitions(){};
  void Init(const G4Fragment & aFragment,G4double dLvlDensity);

  ~PreCompoundTransitions() {}

private:
  
  PreCompoundTransitions(const PreCompoundTransitions &right) {};
  
  const PreCompoundTransitions& operator=(const PreCompoundTransitions &right);

  G4bool operator==(const PreCompoundTransitions &right) const;
  
  G4bool operator!=(const PreCompoundTransitions &right) const;

public:
  G4double GetTotalProbability(void)
  { return TransitionProb1+TransitionProb2+TransitionProb3; }
  
  void PerformTransition(G4Fragment & aFragment,bool bForce);
  
private:
	
  G4double TransitionProb1;
  G4double TransitionProb2;
  G4double TransitionProb3;

};

#endif
