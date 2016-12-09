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
// 080417 Add IsZAApplicable method (return false) by T. Koi
// 080428 Add bool onFlightDB by T. Koi
// 091118 Add Ignore and Enable On Flight Doppler Broadening methods by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPInelasticData_h
#define G4ParticleHPInelasticData_h 1

// Class Description
// Cross-section data set for a high precision (based on evaluated data
// libraries) description of neutron inelastic scattering below 20 MeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"
#include "G4Neutron.hh"

class G4ParticleHPData;

class G4ParticleHPInelasticData : public G4VCrossSectionDataSet
{
   public:
   
      G4ParticleHPInelasticData(G4ParticleDefinition* projectile = G4Neutron::Neutron());
   
      ~G4ParticleHPInelasticData();

      G4bool IsIsoApplicable( const G4DynamicParticle* , 
                              G4int /*Z*/ , G4int /*A*/ ,
                              const G4Element* /*elm*/ ,
                              const G4Material* /*mat*/ );

      G4double GetIsoCrossSection( const G4DynamicParticle*  , 
                                   G4int /*Z*/ , G4int /*A*/ ,
                                   const G4Isotope* /*iso*/  ,
                                   const G4Element* /*elm*/  ,
                                   const G4Material* /*mat*/ );
   
       //G4bool IsApplicable(const G4DynamicParticle*, const G4Element*); ??

   public:
      //G4bool IsZAApplicable( const G4DynamicParticle* , G4double /*ZZ*/, G4double /*AA*/)
      //{ return false;} 

      G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);

  //void BuildPhysicsTableHP(G4ParticleDefinition* projectile, const char* dataDirVariable); // name it 'HP' to avoid compilation warning because of G4VCrossSectionDataSet's method
      void BuildPhysicsTable( const G4ParticleDefinition& ); 

      void DumpPhysicsTable(const G4ParticleDefinition&);

      void IgnoreOnFlightDopplerBroadening(){ onFlightDB = false; };
      void EnableOnFlightDopplerBroadening(){ onFlightDB = true; };
   
  G4ParticleDefinition* GetProjectile(){return theProjectile;}
 
   public:
      G4int GetVerboseLevel() const;
      void SetVerboseLevel( G4int );
      virtual void CrossSectionDescription(std::ostream&) const;
 
   private:
   
      G4PhysicsTable * theCrossSections;

      G4bool onFlightDB;

  G4ParticleDefinition* theProjectile;

      G4ParticleHPData* theHPData;

      G4bool instanceOfWorker;

      G4double ke_cache;
      G4double xs_cache;
      const G4Element* element_cache;
      const G4Material* material_cache;

  //char envVariable[500];
};

#endif
