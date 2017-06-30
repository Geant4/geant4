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
// $Id: G4FissionProbability.hh 103162 2017-03-20 09:40:58Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//
#ifndef G4FissionProbability_h
#define G4FissionProbability_h 1

#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"

class G4FissionProbability : public G4VEmissionProbability
{
public:
  // Default constructor
  explicit G4FissionProbability();

  virtual ~G4FissionProbability();  

  virtual G4double EmissionProbability(const G4Fragment & fragment, 
                                       G4double MaximalKineticEnergy);

  inline void SetEvaporationLevelDensityParameter(G4VLevelDensityParameter * aLevelDensity)
  { 
    if (ownEvapLDP) delete theEvapLDP;
    theEvapLDP = aLevelDensity;
    ownEvapLDP = false;
  }

  inline void SetFissionLevelDensityParameter(G4VLevelDensityParameter * aLevelDensity)
  { 
    if (ownFissLDP) delete theFissLDP;
    theFissLDP = aLevelDensity;
    ownFissLDP = false;
  }

private:

  // Copy constructor
  G4FissionProbability(const G4FissionProbability &right) = delete;
  const G4FissionProbability & operator=(const G4FissionProbability &right) = delete;
  G4bool operator==(const G4FissionProbability &right) const = delete;
  G4bool operator!=(const G4FissionProbability &right) const = delete;
  
  G4VLevelDensityParameter *theEvapLDP;
  G4VLevelDensityParameter *theFissLDP;
  bool ownEvapLDP;
  bool ownFissLDP;


};


#endif
