// by JPW, working, but to be cleaned up. @@@@
// note that the data are taken from carbon12 only.

#ifndef CCalCarbonNonelastic_h
#define CCalCarbonNonelastic_h

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4Proton.hh"

#include "CCalDataSet.hh"

class CCalCarbonNonelastic : public G4VCrossSectionDataSet
{
   public:
   
   CCalCarbonNonelastic();
   virtual ~CCalCarbonNonelastic() {}
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element*anEle)
   {
     G4bool result = false;
     if(( 6 == anEle->GetZ()) &&
        ( aPart->GetKineticEnergy()<150*MeV) &&
	( G4Proton::Proton() == aPart->GetDefinition() ) ) result = true;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle* aPart, const G4Element*, G4double aTemp)
   {
      G4double result = 0;
      const G4double kineticEnergy = aPart->GetKineticEnergy()/eV;
      result = theData.getCrossSection(kineticEnergy);
      return result;
   }

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "CCalCarbonNonelastic: uses ADL data"<<endl;}
   
   private:
   
   CCalDataSet theData;

};

#endif

