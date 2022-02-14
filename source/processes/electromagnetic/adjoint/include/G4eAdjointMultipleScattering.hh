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
// -----------------------------------------------------------------------------
//
// File name:     G4eAdjointMultipleScattering
//
// Author:        Vladimir Ivanchenko
//
// The class simulates the multiple scattering for e+ and e-
//
//------------------------------------------------------------------------------

#ifndef G4eAdjointMultipleScattering_h
#define G4eAdjointMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

class G4eAdjointMultipleScattering : public G4VMultipleScattering

{
 public:
  explicit G4eAdjointMultipleScattering(const G4String& processName = "msc");

  ~G4eAdjointMultipleScattering() override;

  void StartTracking(G4Track*) override;

  // returns true for charged particles, false otherwise
  G4bool IsApplicable(const G4ParticleDefinition& p) override;

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override { ProcessDescription(G4cout); };
  void StreamProcessInfo(std::ostream& outFile) const override;

  G4eAdjointMultipleScattering(G4eAdjointMultipleScattering&) = delete;
  G4eAdjointMultipleScattering& operator                      =(
    const G4eAdjointMultipleScattering& right) = delete;

 protected:
  void InitialiseProcess(const G4ParticleDefinition*) override;

 private:
  G4bool fIsInitialized = false;
};

#endif
