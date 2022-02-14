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
////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointhMultipleScattering
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  The class simulates the multiple scattering for adjoint proton of charged
//  particle. In this approximate implementation the reverse multiple scattering
//  is the same as the forward one. This should be changed in the future to
//  have the MultipleScattering computed for the energy  at the end of the step
//  and not before the step.
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointhMultipleScattering_h
#define G4AdjointhMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

class G4VMscModel;

class G4AdjointhMultipleScattering : public G4VMultipleScattering
{
 public:
  explicit G4AdjointhMultipleScattering(const G4String& processName = "msc");

  ~G4AdjointhMultipleScattering() override;

  // returns true for charged particles, false otherwise
  G4bool IsApplicable(const G4ParticleDefinition& p) override;

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override { ProcessDescription(G4cout); };
  void StreamProcessInfo(std::ostream& out) const override;

  G4AdjointhMultipleScattering(G4AdjointhMultipleScattering&) = delete;
  G4AdjointhMultipleScattering& operator=(
    const G4AdjointhMultipleScattering& right) = delete;

 protected:
  void InitialiseProcess(const G4ParticleDefinition*) override;

 private:
  G4bool fIsInitialized = false;
};

#endif
