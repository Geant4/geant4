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
#ifndef G4ParticleHPThermalScatteringData_h
#define G4ParticleHPThermalScatteringData_h 1

// Thermal Neutron Scattering
// Koi, Tatsumi (SCCS/SLAC)
//
// Class Description
// Cross Sections for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV;
// Based on Thermal neutron scattering files
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

// 15-Nov-06 First implementation is done by T. Koi (SLAC/SCCS)
// 070625 create clearCurrentXSData to fix memory leaking by T. Koi
// 080417 Add IsZAApplicable method (return false) by T. Koi

#include "G4ParticleHPThermalScatteringNames.hh"
#include "G4ParticleHPVector.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
//#include "G4PhysicsTable.hh"

#include <map> 
#include <vector> 

class G4ParticleHPThermalScatteringData : public G4VCrossSectionDataSet
{

   public:

      G4ParticleHPThermalScatteringData();
   
      ~G4ParticleHPThermalScatteringData();

      G4bool IsIsoApplicable( const G4DynamicParticle* , 
                              G4int /*Z*/ , G4int /*A*/ ,
                              const G4Element* /*elm*/ ,
                              const G4Material* /*mat*/ );

      G4double GetIsoCrossSection( const G4DynamicParticle*  , 
                                   G4int /*Z*/ , G4int /*A*/ ,
                                   const G4Isotope* /*iso*/  ,
                                   const G4Element* /*elm*/  ,
                                   const G4Material* /*mat*/ );
   
      G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

      //G4bool IsZAApplicable( const G4DynamicParticle* , G4double /*ZZ*/, G4double /*AA*/)
      //{ return false;}

      G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, const G4Material* );
      G4double GetInelasticCrossSection(const G4DynamicParticle*, const G4Element*, const G4Material*);
      G4double GetCoherentCrossSection(const G4DynamicParticle*, const G4Element*, const G4Material*);
      G4double GetIncoherentCrossSection(const G4DynamicParticle*, const G4Element*, const G4Material*);

      void BuildPhysicsTable(const G4ParticleDefinition&);

      void DumpPhysicsTable(const G4ParticleDefinition&);

      //For user prepared thermal files 
                              //Name of G4Element , Name of NDL file
      void AddUserThermalScatteringFile( G4String , G4String );
   
      virtual void CrossSectionDescription(std::ostream&) const;

   private:

      G4double GetX ( const G4DynamicParticle* , G4double aT , std::map< G4double , G4ParticleHPVector* >* );

      G4double emax; 
   
      void clearCurrentXSData();

//              element            temp       x section from E
      std::map< G4int , std::map< G4double , G4ParticleHPVector* >* >* coherent;
      std::map< G4int , std::map< G4double , G4ParticleHPVector* >* >* incoherent;
      std::map< G4int , std::map< G4double , G4ParticleHPVector* >* >* inelastic;

      std::map< G4double , G4ParticleHPVector* >* readData ( G4String ); 

      std::vector < G4int > indexOfThermalElement;
      G4ParticleHPThermalScatteringNames* names;
//              G4Element  NDL name 
//      std::map< G4String , G4String > names;

      G4double ke_cache;
      G4double xs_cache;
      const G4Element* element_cache;
      const G4Material* material_cache;

      std::map < std::pair < const G4Material* , const G4Element* > , G4int > dic;   
      G4int getTS_ID( const G4Material* , const G4Element* );

};

#endif
