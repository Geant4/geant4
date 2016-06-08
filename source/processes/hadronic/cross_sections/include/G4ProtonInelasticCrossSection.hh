// by JPW, working, but to be cleaned up. @@@@

#ifndef G4ProtonInelasticCrossSection_h
#define G4ProtonInelasticCrossSection_h

// Class Description
// Cross-sections for proton nuclear scattering up to 20 GeV, getting the low
// energy threshold behaviour right.
// H.P. Wellisch (TRIUMF), D. Axen (British Columbia U.). 1996. 
// Published in Phys.Rev.C54:1329-1332,1996 
// Class Description - End
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
   
   G4double GetCrossSection(G4double anEnergy, G4double anA, G4double aZ);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4ProtonInelasticCrossSection: uses formula"<<G4endl;}

};

#endif
