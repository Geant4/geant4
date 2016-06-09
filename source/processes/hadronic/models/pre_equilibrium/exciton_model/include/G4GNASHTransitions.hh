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
#ifndef G4GNASHTransitions_hh
#define G4GNASHTransitions_hh 1

#include "G4VPreCompoundTransitions.hh"
#include "G4Fragment.hh"

class G4GNASHTransitions : public G4VPreCompoundTransitions
{
public:
  G4GNASHTransitions() {}
  virtual ~G4GNASHTransitions() {}

private:

  G4GNASHTransitions(const G4GNASHTransitions &) : G4VPreCompoundTransitions() {}

  const G4GNASHTransitions& operator=(const G4GNASHTransitions &right);

  G4bool operator==(const G4GNASHTransitions &right) const;

  G4bool operator!=(const G4GNASHTransitions &right) const;

public:

  virtual G4double CalculateProbability(const G4Fragment & aFragment);
  
  virtual G4Fragment PerformTransition(const G4Fragment & aFragment);


};

#endif
