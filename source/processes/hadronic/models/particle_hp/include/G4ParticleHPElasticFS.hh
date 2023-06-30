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
//
// 25-08-06 New Final State type (refFlag==3 , Legendre (Low Energy) + Probability (High Energy) )
//          is added by T. KOI
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPElasticFS_h
#define G4ParticleHPElasticFS_h 1

#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4ParticleHPFastLegendre.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPInterpolator.hh"
#include "G4ParticleHPLegendreStore.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPPartial.hh"
#include "G4ParticleHPVector.hh"
#include "G4ReactionProduct.hh"
#include "globals.hh"

class G4ParticleHPElasticFS : public G4ParticleHPFinalState
{
  public:
    G4ParticleHPElasticFS();
    ~G4ParticleHPElasticFS() override
    {
      delete theCoefficients;
      delete theProbArray;
    }
    void Init(G4double A, G4double Z, G4int M, G4String& dirName, G4String& aFSType,
              G4ParticleDefinition*) override;
    G4HadFinalState* ApplyYourself(const G4HadProjectile& theTrack) override;
    G4ParticleHPFinalState* New() override
    {
      auto theNew = new G4ParticleHPElasticFS;
      return theNew;
    }

    // New method useful for the DBRC algorithm
    void RegisterCrossSection(G4ParticleHPVector* vec) { xsForDBRC = vec; }

  private:
    // The following two methods and six data members are needed for the DBRC algorithm
    void InitializeScatteringKernelParameters();
    G4ReactionProduct GetBiasedThermalNucleus(const G4double aMass, G4ThreeVector aVelocity,
                                              const G4double temp = -1.);
    G4double svtEmax;
    G4double dbrcEmax;
    G4double dbrcEmin;
    G4double dbrcAmin;
    G4bool dbrcUse;
    G4ParticleHPVector* xsForDBRC;

    G4int repFlag;  // Legendre coeff(1), or probability array(2), or isotropic(0).
                    // add 3: for Legendre (Low Energy) + Probability (High Energy)
    G4double tE_of_repFlag3;  // transition energy in case of  repFlag 3:
    G4double targetMass;  // in neutronmass units.
    G4int frameFlag;  // CMS or Lab system.

    G4ParticleHPLegendreStore* theCoefficients;  // the legendre coefficients
    G4ParticleHPPartial* theProbArray;  // the probability array p,costh for energy

    G4ParticleHPFastLegendre theLegend;  // fast look-up for leg-integrals
};
#endif
