// by JPW, working, but to be cleaned up. @@@@

#ifndef G4TripathiCrossSection_h
#define G4TripathiCrossSection_h

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4TripathiCrossSection : public G4VCrossSectionDataSet
{
   public:
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element*)
   {
     G4bool result = false;
     if(( aPart->GetDefinition()->GetBaryonNumber()>2.5) &&
        ( aPart->GetKineticEnergy()/aPart->GetDefinition()->GetBaryonNumber()<1*GeV) ) result = true;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4TripathiCrossSection: uses formula"<<G4endl;}

};

#endif
