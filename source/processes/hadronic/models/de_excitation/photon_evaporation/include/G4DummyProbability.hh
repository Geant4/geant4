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
//

#ifndef G4DummyProbability_hh
#define G4DummyProbability_hh

#include "globals.hh"
#include "G4VEmissionProbability.hh"
#include "G4Fragment.hh"
#include "G4VLevelDensityParameter.hh"

class G4DummyProbability : public G4VEmissionProbability
{

public:

  G4DummyProbability() {};

  ~G4DummyProbability();

  G4double EmissionProbability(const G4Fragment& frag, const G4double excite);
  G4double EmissionProbDensity(const G4Fragment& frag, const G4double ePhoton);

private:

  // G4DummyProbability() {};

  G4DummyProbability(const G4DummyProbability& right);

  const G4DummyProbability& operator=(const G4DummyProbability& right);
  G4bool operator==(const G4DummyProbability& right) const;
  G4bool operator!=(const G4DummyProbability& right) const;

};

#endif





