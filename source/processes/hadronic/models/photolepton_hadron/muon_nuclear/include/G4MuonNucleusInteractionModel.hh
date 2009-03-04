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
// G4MuonNucleusInteractionModel.hh
//
//     M.Takahata (Makoto.Takahata@cern.ch)
//
// -----------------------------------------------------------------------
//  This class is the translation of Geant3/Gheisha routine GMUNU + GMUSIG
//  , which calculate the final state of outgoing muon and force cascade
//  through pion-Nucleus inelastic scattering.
// -----------------------------------------------------------------------

#ifndef G4MuonNucleusInteractionModel_h
#define G4MuonNucleusInteractionModel_h 1

#include "G4LeptonHadronInteractionModel.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicsLogVector.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HadFinalState.hh"


  class G4MuonNucleusInteractionModel : public G4LeptonHadronInteractionModel
  {
    public:

      G4MuonNucleusInteractionModel();
      virtual ~G4MuonNucleusInteractionModel();

      void makePhysicsVector();
      G4double computeMicroscopicCrossSection(const G4Track &muonTrack);
      G4VParticleChange* applyInteractionModel(const G4Track &muonTrack, 
                                               G4Nucleus &targetNucleus);

    private:

      // member functions
      void invokePionNucleus(const G4Track &pionTrack,
                             G4Nucleus &targetNucleus);
      G4double computeDifferentialCrossSection
        (G4double initialEnergy, G4double finalEnergy, G4double cosTheta);

      // for physics vector
      static G4double tetal[35], xeml[23];
      G4PhysicsLogVector *theCoefficientVector;
      G4int    Nbin;
      G4double kEmin, kEmax;

      // kinematic
      G4double cosTheta;
      G4HadFinalState  * pionChange;

      // pi-N models
      G4double cascadeModelMarginalEnergy;
      G4LEPionPlusInelastic*  LEPionPlusInelastic;
      G4LEPionMinusInelastic* LEPionMinusInelastic;
      G4HEPionPlusInelastic*  HEPionPlusInelastic;
      G4HEPionMinusInelastic* HEPionMinusInelastic;
  };

#endif
