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
// $Id: G4AngleDirect.hh 73844 2013-09-13 14:16:30Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4AngleDirect
//
// Author:        V. Ivanchenko 
// 
// Creation date: 03 October 2013
//
// Modifications: 
//
// Class Description: 
//
// Secondary particle has the same direction as a primary
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4AngleDirect_h
#define G4AngleDirect_h 1

#include "G4VEmAngularDistribution.hh"

class G4Material;

class G4AngleDirect : public G4VEmAngularDistribution
{
public:

  G4AngleDirect();

  virtual ~G4AngleDirect();

  // Sample direction in global coordinate system,
  // this means for zero scattering angle this direction is the same
  // as the direction of primary 
  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
					 G4double, G4int, 
					 const G4Material*) override;

private:

  // hide assignment operator 
  G4AngleDirect & operator=(const  G4AngleDirect &right) = delete;
  G4AngleDirect(const  G4AngleDirect&) = delete;

};

#endif

