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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4VSubCutProducer
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

#ifndef G4VSubCutProducer_h
#define G4VSubCutProducer_h 1

#include "globals.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include <vector>

class G4VSubCutProducer
{
public:

  explicit G4VSubCutProducer(const G4String& name) : fName(name) {};

  virtual ~G4VSubCutProducer() = default;

  // Sample direction in global coordinate system,
  // this means for zero scattering angle this direction is the same
  // as the direction of primary 
  virtual void SampleSecondaries(const G4Step& step,
				 std::vector<G4Track*>& tracks,
                                 G4double& eloss,
				 G4double cut) const = 0;

  inline const G4String& GetName() const;

  // hide assignment operator 
  G4VSubCutProducer & operator=(const  G4VSubCutProducer &right) = delete;
  G4VSubCutProducer(const  G4VSubCutProducer&) = delete;

private:

  G4String fName;
};

inline const G4String& G4VSubCutProducer::GetName() const
{
  return fName;
}

#endif

