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
// $Id: G4AtomicTransition.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// ?????      Created
//
// -------------------------------------------------------------------

#ifndef G4AtomicTransition_h 
#define G4AtomicTransition_h 1

#include "G4DataVector.hh"
#include "globals.hh"
#include "g4std/vector"

class G4AtomicTransition {

public:

  G4AtomicTransition(G4int,const G4std::vector<G4int>&,const G4DataVector&,
		     const G4DataVector&);
  ~G4AtomicTransition();

  const G4std::vector<G4int>& OriginatingShellIds() const;
  const G4DataVector& TransitionEnergies() const;
  const G4DataVector& TransitionProbabilities() const;
  const G4int FinalShellId() const;
  G4int OriginatingShellId(G4int index) const;
  G4double TransitionEnergy(G4int index) const;
  G4double TransitionProbability(G4int index) const;

private:

  G4int finalShellId;
  G4std::vector<G4int> originatingShellIds;
  G4DataVector transitionEnergies;
  G4DataVector transitionProbabilities;
  
};

#endif

