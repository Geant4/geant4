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
// $Id: G4PreCompoundTransitions.hh 89523 2015-04-16 09:56:35Z gcosmo $
//
// by V. Lara
// 01.05.2008 J. M. Quesada . New methods for accessing to individual transition 
//                 probabilities (landa+, landa-, landa0) from G4PreCompoundModel  
// 20.08.2010 V.Ivanchenko move constructor and destructor to the source

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

class G4ParticleDefinition;
class G4Fragment;

class G4PreCompoundTransitions : public G4VPreCompoundTransitions
{
public:

  G4PreCompoundTransitions();

  virtual ~G4PreCompoundTransitions();

  virtual G4double CalculateProbability(const G4Fragment & aFragment);
  
  virtual void PerformTransition(G4Fragment & aFragment);

private:
  
  G4PreCompoundTransitions(const G4PreCompoundTransitions &); 
  const G4PreCompoundTransitions& operator=(const G4PreCompoundTransitions &right);
  G4bool operator==(const G4PreCompoundTransitions &right) const;
  G4bool operator!=(const G4PreCompoundTransitions &right) const;

  const G4ParticleDefinition* proton;

  G4double FermiEnergy;
  G4double r0;  // Nuclear radius
  G4double aLDP;// Level density parameter
};

#endif
