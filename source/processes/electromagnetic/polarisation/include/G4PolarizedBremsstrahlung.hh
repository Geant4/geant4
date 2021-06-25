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
// -------------------------------------------------------------------
//
// Geant4 Class header file
//
// File name:     G4PolarizedBremsstrahlung
//
// Author:        Karim Laihem based on code by Michel Maire
//
// Class Description:
//   polarized version of G4eBremsstrahlung

#ifndef G4PolarizedBremsstrahlung_h
#define G4PolarizedBremsstrahlung_h 1

#include "globals.hh"
#include "G4eBremsstrahlung.hh"

class G4ParticleDefinition;

class G4PolarizedBremsstrahlung : public G4eBremsstrahlung
{
 public:
  explicit G4PolarizedBremsstrahlung(const G4String& name = "pol-eBrem");
  virtual ~G4PolarizedBremsstrahlung() override;

  virtual void ProcessDescription(std::ostream&) const override;
  virtual void DumpInfo() const override { ProcessDescription(G4cout); };

  G4PolarizedBremsstrahlung& operator=(const G4PolarizedBremsstrahlung& right) =
    delete;
  G4PolarizedBremsstrahlung(const G4PolarizedBremsstrahlung&) = delete;

 protected:
  virtual void InitialiseEnergyLossProcess(
    const G4ParticleDefinition*, const G4ParticleDefinition*) override;

private:
  G4bool isInitialised = false;
};
#endif
