// by JPW, working, but to be cleaned up. @@@@

#ifndef G4ProtonInelasticCrossSection_h
#define G4ProtonInelasticCrossSection_h

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4ProtonInelasticCrossSection : public G4VCrossSectionDataSet
{
   public:
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element*)
   {
     G4bool result = false;
     if(( aPart->GetDefinition()==G4Proton::Proton()) &&
        ( aPart->GetKineticEnergy()<20*GeV) ) result = true;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4ProtonInelasticCrossSection: uses formula"<<G4endl;}

};

#endif
