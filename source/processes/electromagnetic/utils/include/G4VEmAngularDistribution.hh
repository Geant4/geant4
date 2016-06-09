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
// $Id: G4VEmAngularDistribution.hh,v 1.2 2010-11-04 12:55:09 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4VEmAngularDistribution
{
public:

  G4VEmAngularDistribution(const G4String& name);

  virtual ~G4VEmAngularDistribution();

  // method for bremsstrahlung
  virtual G4double PolarAngle(const G4double initial_energy,
			      const G4double final_energy,
                              const G4int Z) = 0;

  inline const G4String& GetName() const;

private:

  // hide assignment operator 
  G4VEmAngularDistribution & operator=(const  G4VEmAngularDistribution &right);
  G4VEmAngularDistribution(const  G4VEmAngularDistribution&);

  G4String fName;
};

inline const G4String& G4VEmAngularDistribution::GetName() const
{
  return fName;
}

#endif

