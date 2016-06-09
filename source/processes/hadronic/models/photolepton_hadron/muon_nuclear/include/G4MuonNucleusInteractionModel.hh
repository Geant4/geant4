//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
      ~G4MuonNucleusInteractionModel();

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
