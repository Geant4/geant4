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
// $Id: G4Vee2hadrons.hh,v 1.4 2008-07-10 18:06:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  G4Vee2hadrons() : lowEnergy(0.0), highEnergy(1.1*GeV) {};

  virtual ~G4Vee2hadrons() {};

  virtual G4double ThresholdEnergy() const = 0;

  virtual G4double PeakEnergy() const = 0;

  virtual G4double ComputeCrossSection(G4double) const = 0;

  virtual G4PhysicsVector* PhysicsVector(G4double, G4double) const = 0;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 G4double, const G4ThreeVector&) = 0;

  void SetLowEnergy(G4double val) {lowEnergy = val;};

  G4double LowEnergy() const {return lowEnergy;};

  void SetHighEnergy(G4double val) {highEnergy = val;};

  G4double HighEnergy() const {return highEnergy;};

private:

  // hide assignment operator
  G4Vee2hadrons & operator=(const  G4Vee2hadrons &right);
  G4Vee2hadrons(const  G4Vee2hadrons&);

  G4double lowEnergy;
  G4double highEnergy;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
