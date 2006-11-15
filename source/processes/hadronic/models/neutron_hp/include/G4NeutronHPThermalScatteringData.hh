#ifndef G4NeutronHPThermalScatteringData_h
#define G4NeutronHPThermalScatteringData_h 1

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

#include "G4NeutronHPThermalScatteringNames.hh"
#include "G4NeutronHPVector.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
//#include "G4PhysicsTable.hh"

#include <map> 
#include <vector> 

class G4NeutronHPThermalScatteringData : public G4VCrossSectionDataSet
{
   public:
   
      G4NeutronHPThermalScatteringData();
   
      ~G4NeutronHPThermalScatteringData();
   
      G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

      G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);
      G4double GetInelasticCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);
      G4double GetCoherentCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);
      G4double GetIncoherentCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);

      void BuildPhysicsTable(const G4ParticleDefinition&);

      void DumpPhysicsTable(const G4ParticleDefinition&);
   
   private:

      G4double GetX ( const G4DynamicParticle* , G4double aT , std::map< G4double , G4NeutronHPVector* >* );

      G4double emax; 
   

//              element            temp       x section from E
      std::map< G4int , std::map< G4double , G4NeutronHPVector* >* > coherent;
      std::map< G4int , std::map< G4double , G4NeutronHPVector* >* > incoherent;
      std::map< G4int , std::map< G4double , G4NeutronHPVector* >* > inelastic;

      std::map< G4double , G4NeutronHPVector* >* readData ( G4String ); 

      std::vector < G4int > indexOfThermalElement;
      G4NeutronHPThermalScatteringNames* names;
//              G4Element  NDL name 
//      std::map< G4String , G4String > names;

//   G4PhysicsTable * theCrossSections;

};

#endif
