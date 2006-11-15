#ifndef G4NeutronHPJENDLHEData_h
#define G4NeutronHPJENDLHEData_h 1

// Class Description
// Cross-section data set for a high precision (based on JENDL_HE evaluated data
// libraries) description of elastic scattering 20 MeV ~ 3 GeV; 
// Class Description - End

// 15-Nov-06 First Implementation is done by T. Koi (SLAC/SCCS)

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include <map> 

class G4NeutronHPJENDLHEData : public G4VCrossSectionDataSet
{
   public:
   
   G4NeutronHPJENDLHEData();
   G4NeutronHPJENDLHEData( G4String , G4ParticleDefinition* );

   
   ~G4NeutronHPJENDLHEData();
   
   G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);

   void BuildPhysicsTable(const G4ParticleDefinition&);

   void DumpPhysicsTable(const G4ParticleDefinition&);
   
   private:
   
      std::vector< G4bool > vElement;

      std::map< G4int , std::map< G4int , G4PhysicsVector* >* > mIsotope; 

      G4bool isThisInMap ( G4int , G4int );
      G4bool isThisNewIsotope ( G4int z , G4int a ) { return !( isThisInMap( z , a ) ); };
      G4PhysicsVector* readAFile ( std::fstream* );
      void registAPhysicsVector ( G4int , G4int , G4PhysicsVector* );

      G4double getXSfromThisIsotope ( G4int , G4int , G4double );

      G4String reactionName;
      G4String particleName;
};

#endif
