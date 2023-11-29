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
// The model of transition radiation
//
// History:
//
// 04.10.05, V.Grichine move from pure virtual and new class name
// 29.02.04, V.Ivanchenko created

#ifndef G4VTRModel_h
#define G4VTRModel_h

#include "globals.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"

#include <vector>

class G4Track;

class G4VTRModel
{
 public:
  // Constructors
  explicit G4VTRModel(const G4String& modelName) { fName = modelName; };

  // Destructor
  virtual ~G4VTRModel() = default;

  const G4String& GetName() const { return fName; };

  virtual void GenerateSecondaries(G4VParticleChange& pChange,
                                   std::vector<const G4Material*>& materials,
                                   std::vector<G4double>& steps,
                                   std::vector<G4ThreeVector>& normals,
                                   G4ThreeVector& startingPosition,
                                   const G4Track& track);

  // disable assignment operator & copy constructor
  G4VTRModel& operator=(const G4VTRModel& right) = delete;
  G4VTRModel(const G4VTRModel&)                  = delete;

  virtual void PrintInfo() { return; };

 protected:
  G4String fName;
};

#endif  // G4VTRModel_h
