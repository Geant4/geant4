// by JPW, working, but to be cleaned up. @@@@

#ifndef G4NeutronInelasticCrossSection_h
#define G4NeutronInelasticCrossSection_h

#include "globals.hh"
#include "G4Neutron.hh"

#include "G4VCrossSectionDataSet.hh"

class G4NeutronInelasticCrossSection : public G4VCrossSectionDataSet
{
   public:
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element*)
   {
     G4bool result = false;
     if(( aPart->GetDefinition()==G4Neutron::Neutron()) &&
        ( aPart->GetKineticEnergy()<20*GeV) &&
          aPart->GetKineticEnergy()>19.9*MeV) result = true;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4NeutronInelasticCrossSection: uses formula"<<endl;}

};

#endif
