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
#ifndef G4ParticleHPThermalScattering_h
#define G4ParticleHPThermalScattering_h 1

// Thermal Neutron Scattering
// Koi, Tatsumi (SLAC/SCCS)
//
// Class Description
// Final State Generators for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV;
// Based on Thermal neutron scattering files
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

#include "globals.hh"
#include "G4ParticleHPThermalScatteringNames.hh"
#include "G4HadronicInteraction.hh"

class G4ParticleHPThermalScatteringData;
class G4ParticleHPElastic;

struct E_isoAng 
{
   G4double energy;
   G4int n;
   std::vector < G4double > isoAngle; 
   E_isoAng() {
      energy=0.0;
      n=0;
   };
};

struct E_P_E_isoAng 
{
   G4double energy;
   G4int n;
   std::vector < G4double > prob; 
   std::vector < E_isoAng* > vE_isoAngle; 
   G4double sum_of_probXdEs;  // should be close to 1
   std::vector< G4double > secondary_energy_cdf;
   std::vector< G4double > secondary_energy_pdf;
   std::vector< G4double > secondary_energy_value; 
   G4int secondary_energy_cdf_size;
   E_P_E_isoAng() {
      energy=0.0;
      n=0;
      sum_of_probXdEs=0.0;
      secondary_energy_cdf_size=0;
   };
};

class G4ParticleHPThermalScattering : public G4HadronicInteraction
{
   public: 
  
      G4ParticleHPThermalScattering();
  
      ~G4ParticleHPThermalScattering();
  
      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);

      virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

      //For user prepared thermal files 
                              //Name of G4Element , Name of NDL file
      void AddUserThermalScatteringFile( G4String , G4String );

      void BuildPhysicsTable(const G4ParticleDefinition&);

      virtual void ModelDescription(std::ostream& outFile) const;

   private:

      void clearCurrentFSData();

      G4ParticleHPThermalScatteringNames names;

      // Coherent Elastic 
      //         ElementID             temp                             BraggE     cumulativeP
      std::map < G4int , std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* >* coherentFSs;
      std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* readACoherentFSDATA( G4String );

      // Incoherent Elastic 
      //         ElementID          temp       aFS for this temp (and this element)
      std::map < G4int , std::map < G4double , std::vector < E_isoAng* >* >* >* incoherentFSs;
      std::map < G4double , std::vector < E_isoAng* >* >* readAnIncoherentFSDATA( G4String );
      E_isoAng* readAnE_isoAng ( std::istream* );

      // Inelastic 
      //         ElementID          temp        aFS for this temp (and this element) 
      std::map < G4int ,  std::map < G4double , std::vector < E_P_E_isoAng* >* >* >* inelasticFSs;
      std::map < G4double , std::vector < E_P_E_isoAng* >* >* readAnInelasticFSDATA( G4String );
      E_P_E_isoAng* readAnE_P_E_isoAng ( std::istream* );
  
      G4ParticleHPThermalScatteringData* theXSection;

      G4ParticleHPElastic* theHPElastic;
      
      G4double getMu ( E_isoAng* );
      G4double getMu ( G4double rndm1 , G4double rndm2 , E_isoAng* anEPM );

      std::pair< G4double , G4double > find_LH ( G4double , std::vector<G4double>* );
      G4double get_linear_interpolated ( G4double , std::pair < G4double , G4double > , std::pair < G4double , G4double > );

      E_isoAng create_E_isoAng_from_energy( G4double , std::vector< E_isoAng* >* );

      G4double get_secondary_energy_from_E_P_E_isoAng ( G4double random , E_P_E_isoAng* anE_P_E_isoAng );

      std::pair< G4double, G4double > sample_inelastic_E_mu( G4double pE , std::vector< E_P_E_isoAng* >* vNEP_EPM ); 
      std::pair< G4double, G4int > sample_inelastic_E( G4double rndm1 , G4double rndm2 , E_P_E_isoAng* anE_P_E_isoAng );

      std::pair< G4double , E_isoAng > create_sE_and_EPM_from_pE_and_vE_P_E_isoAng ( G4double , G4double , std::vector < E_P_E_isoAng* >* );

      std::map < std::pair < const G4Material* , const G4Element* > , G4int > dic;   
      void buildPhysicsTable();
      G4int getTS_ID( const G4Material* , const G4Element* );

      //size_t sizeOfMaterialTable;

      G4bool check_E_isoAng( E_isoAng* );

      //In order to judge whether the rebuilding of physics table is a necessity or not
      size_t nMaterial;
      size_t nElement;

};

#endif
