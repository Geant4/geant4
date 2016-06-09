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
// $Id: G4Vee2hadrons.hh,v 1.1 2004/11/19 18:44:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4Vee2hadrons
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 12.08.2004
//
// Modifications:
//

//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4Vee2hadrons_h
#define G4Vee2hadrons_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4DynamicParticle;
class G4PhysicsVector;

class G4Vee2hadrons 
{

public:

  G4Vee2hadrons() {};

  virtual ~G4Vee2hadrons() {};

  virtual G4double ThresholdEnergy() const = 0;

  virtual G4double PeakEnergy() const = 0;

  virtual G4double ComputeCrossSection(G4double) const = 0;

  virtual G4PhysicsVector* PhysicsVector(G4double, G4double) const = 0;

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                      G4double, const G4ThreeVector&) const = 0;

private:

  // hide assignment operator
  G4Vee2hadrons & operator=(const  G4Vee2hadrons &right);
  G4Vee2hadrons(const  G4Vee2hadrons&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
