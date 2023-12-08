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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPIsoData_h
#define G4ParticleHPIsoData_h 1

// Hadronic Process: Very Low Energy Neutron X-Sections
// original by H.P. Wellisch, TRIUMF, 14-Feb-97
// Has the Cross-section data for on isotope.

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPVector.hh"

class G4ParticleDefinition;

class G4ParticleHPIsoData
{
  public:
    G4ParticleHPIsoData() = default;

    ~G4ParticleHPIsoData() { delete theChannelData; }

    inline G4double GetXsec(G4double energy)
    {
      return std::max(0., theChannelData->GetXsec(energy));
    }

    G4bool Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType)
    {
      return Init(A, Z, 0, abun, dirName, aFSType);
    };

    G4bool Init(G4int A, G4int Z, G4int M, G4double abun, G4String dirName, G4String aFSType);

    void Init(G4int A, G4int Z, G4double abun, G4ParticleDefinition* projectile,
              const char* dataDirVariable)
    {
      Init(A, Z, 0, abun, projectile, dataDirVariable);
    }

    void Init(G4int A, G4int Z, G4int M, G4double abun, G4ParticleDefinition* projectile,
              const char* dataDirVariable);  // fill PhysicsVector for this Isotope

    G4ParticleHPVector* MakeElasticData() { return theElasticData; }
    G4ParticleHPVector* MakeFissionData() { return theFissionData; }
    G4ParticleHPVector* MakeCaptureData() { return theCaptureData; }
    G4ParticleHPVector* MakeInelasticData() { return theInelasticData; }
    G4ParticleHPVector* MakeChannelData() { return theChannelData; }

    G4String GetName(G4int A, G4int Z, G4String base, G4String rest);

    void FillChannelData(G4ParticleHPVector* aBuffer);

    void ThinOut(G4double precision)
    {
      if (theFissionData != nullptr) theFissionData->ThinOut(precision);
      if (theCaptureData != nullptr) theCaptureData->ThinOut(precision);
      if (theElasticData != nullptr) theElasticData->ThinOut(precision);
      if (theInelasticData != nullptr) theInelasticData->ThinOut(precision);
    }

    G4ParticleHPIsoData& operator=(const G4ParticleHPIsoData& right) = delete;
    G4ParticleHPIsoData(const G4ParticleHPIsoData&) = delete;

  private:
    G4ParticleHPVector* theFissionData{nullptr};
    G4ParticleHPVector* theCaptureData{nullptr};
    G4ParticleHPVector* theElasticData{nullptr};
    G4ParticleHPVector* theInelasticData{nullptr};
    G4ParticleHPVector* theChannelData{nullptr};

    G4String theFileName;
    G4ParticleHPNames theNames;
};

#endif
