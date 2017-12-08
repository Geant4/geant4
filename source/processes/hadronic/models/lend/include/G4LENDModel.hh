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
#ifndef G4LENDModel_h
#define G4LENDModel_h 1

// Class Description
// Final state production model for a LEND (Low Energy Nuclear Data) 
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// LEND is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDHeader.hh"
#include "G4LENDManager.hh"
#include "G4LENDUsedTarget.hh"
#include "G4HadronicInteraction.hh"
#include "globals.hh"

extern "C" double MyRNG(void*); 

class G4LENDModel : public G4HadronicInteraction
{

   public: 
  
      //G4LENDModel();
      G4LENDModel( G4String name="LENDModel" );
      ~G4LENDModel();
  
      virtual G4HadFinalState * ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );

      void ChangeDefaultEvaluation( G4String name ){ default_evaluation = name; recreate_used_target_map(); };
      void AllowNaturalAbundanceTarget(){ allow_nat = true; recreate_used_target_map(); };
      void AllowAnyCandidateTarget(){ allow_any = true; recreate_used_target_map(); };
      //Same argument to the CrossSectionDataSet 
      void BuildPhysicsTable( const G4ParticleDefinition& ){ recreate_used_target_map(); };
      void DumpLENDTargetInfo( G4bool force = false );

   private:

      G4String default_evaluation;
      G4bool allow_nat;
      G4bool allow_any;

   protected:

      void create_used_target_map();
      void recreate_used_target_map();
      G4GIDI_target* get_target_from_map( G4int nuclear_code );

      G4HadFinalState* returnUnchanged( const G4HadProjectile& aTrack, G4HadFinalState* theResult );

      G4ParticleDefinition* proj;
      G4LENDManager* lend_manager;
      std::map< G4int , G4LENDUsedTarget* > usedTarget_map;
};

#endif
