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
// $Id: G4VEmAngularDistribution.hh 95657 2016-02-17 13:03:36Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4VEmAngularDistribution
//
// Author:        V. Ivanchenko using design of existing 
//                interface G4VBremAngularDistribution
// 
// Creation date: 13 October 2010
//
// Modifications: 
//
// Class Description: 
//
// Abstract base class for polar angle sampling
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4VEmAngularDistribution_h
#define G4VEmAngularDistribution_h 1

#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"

class G4Material;

class G4VEmAngularDistribution
{
public:

  explicit G4VEmAngularDistribution(const G4String& name);

  virtual ~G4VEmAngularDistribution();

  // Sample direction in global coordinate system,
  // this means for zero scattering angle this direction is the same
  // as the direction of primary 
  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
					 G4double finalTotalEnergy,
					 G4int Z,
					 const G4Material*) = 0;

  // Sample direction in global coordinate system for given electronic shell,
  // this means for zero scattering angle this direction is the same
  // as the direction of primary 
  virtual G4ThreeVector& SampleDirectionForShell(
					 const G4DynamicParticle* dp,
					 G4double finalTotalEnergy,
					 G4int Z,
					 G4int shellID,
					 const G4Material*);

  inline const G4String& GetName() const;

protected:

  G4ThreeVector fLocalDirection;

private:

  // hide assignment operator
  G4VEmAngularDistribution & 
    operator=(const  G4VEmAngularDistribution &right) = delete;
  G4VEmAngularDistribution(const  G4VEmAngularDistribution&) = delete;

  G4String fName;
};

inline const G4String& G4VEmAngularDistribution::GetName() const
{
  return fName;
}

#endif

