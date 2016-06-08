// by JPW, working, but to be cleaned up. @@@@

#ifndef G4NeutronInelasticCrossSection_h
#define G4NeutronInelasticCrossSection_h

#include "globals.hh"
#include "G4Neutron.hh"

// Class Description
// Cross-sections for neutron nuclear scattering from 14 MeV up to 20 GeV, getting the
// low energy threshold behaviour right.
// H.P. Wellisch (TRIUMF), M. Laidlaw (British Columbia U.). 1996. 
// Class Description - End

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
   
   G4double GetCrossSection(G4double anEnergy, G4double anA, G4double aZ);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4NeutronInelasticCrossSection: uses formula"<<G4endl;}

};

#endif
