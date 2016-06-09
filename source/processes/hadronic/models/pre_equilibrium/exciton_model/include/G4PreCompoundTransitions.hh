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
// $Id: G4PreCompoundTransitions.hh,v 1.3 2006/06/29 20:58:36 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
