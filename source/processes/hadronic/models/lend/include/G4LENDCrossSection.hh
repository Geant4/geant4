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
#ifndef G4LENDCrossSection_h
#define G4LENDCrossSection_h 1

// Class Description
// Cross Sections for a LEND (Low Energy Nuclear Data)
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDHeader.hh"

#include "G4LENDManager.hh"
#include "G4LENDUsedTarget.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"

#include <map>

class G4LENDCrossSection : public G4VCrossSectionDataSet
{

//
//
//  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
//                         const G4Element* elm = 0,
//                         const G4Material* mat = 0);
// 
//  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
//                              const G4Element* elm = 0,
//                              const G4Material* mat = 0);
//

   public:
   
      G4LENDCrossSection(const G4String name = "" );
   
      ~G4LENDCrossSection();
   
      //G4bool IsApplicable( const G4DynamicParticle* , const G4Element* );
//TK110811
      G4bool IsIsoApplicable( const G4DynamicParticle* , G4int /*Z*/ , G4int /*A*/ , 
                              const G4Element* , const G4Material* );
      G4double GetIsoCrossSection( const G4DynamicParticle*, G4int /*Z*/, G4int /*A*/ ,
                                   const G4Isotope* , const G4Element* , const G4Material* );

      void BuildPhysicsTable( const G4ParticleDefinition& );

      void DumpPhysicsTable( const G4ParticleDefinition& );
      void DumpLENDTargetInfo( G4bool force = false );

//TK110810
      //G4double GetCrossSection( const G4DynamicParticle* , const G4Element* , G4double aT );
      //G4double GetCrossSection(const G4DynamicParticle*, G4int , const G4Material* );


      void ChangeDefaultEvaluation( G4String name_tmp ){ default_evaluation = name_tmp; };
      void AllowNaturalAbundanceTarget(){ allow_nat = true; };
      void AllowAnyCandidateTarget(){ allow_any = true; };

      //Hadronic Framework still does not handle isotope in GPIL
      //G4VDiscreteProcess::PostStepGetPhysicalInteractionLenght()
      // G4HadronicProcess::GetMeanFreePath()
      //  G4double GetCrossSection(const G4DynamicParticle*, const G4Material*)
      //   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aTemperature); 

//TK110810
      //G4bool IsIsoApplicable( const G4DynamicParticle* , G4int /*ZZ*/, G4int /*AA*/) { return true; }
      //G4bool IsZandAApplicable( const G4DynamicParticle* , G4int /*ZZ*/, G4int /*AA*/, const G4Element* , const G4Material* ) { return true; }
//TK110810
      //G4double GetZandACrossSection(const G4DynamicParticle* , G4int /*Z*/, G4int /*A*/, G4double aTemperature);
      //G4double GetZandACrossSection(const G4DynamicParticle* , G4int /*Z*/, G4int /*A*/,  const G4Material* mat);
//TK110810
      //G4double GetIsoCrossSection( const G4DynamicParticle* , const G4Isotope* , G4double );
//      G4double GetIsoCrossSection( const G4DynamicParticle* , const G4Isotope* , const G4Material* mat );

   private:

      std::map< G4int , G4LENDUsedTarget* > usedTarget_map;

      G4String default_evaluation;
      G4bool allow_nat;
      G4bool allow_any;

      G4LENDManager* lend_manager;
      void recreate_used_target_map();

   protected :
      //G4String name;
      G4ParticleDefinition* proj;
      void create_used_target_map();
      G4GIDI_target* get_target_from_map( G4int nuclear_code );
      virtual G4double getLENDCrossSection( G4GIDI_target* , G4double , G4double ) { return 0.0; };
      //                                         elow       ehigh     xs_elow   xs_ehigh    ke (<elow)
      G4double GetUltraLowEnergyExtrapolatedXS( G4double , G4double , G4double , G4double , G4double );
};
#endif
