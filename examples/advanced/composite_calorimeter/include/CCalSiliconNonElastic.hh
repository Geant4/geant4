#ifndef CCalSiliconNonElastic_h
#define CCalSiliconNonElastic_h

// by JPW, working, but to be cleaned up. @@@@

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4Proton.hh"

#include "CCalDataSet.hh"

class CCalSiliconNonElastic : public G4VCrossSectionDataSet
{
   public:
   
   CCalSiliconNonElastic();
   virtual ~CCalSiliconNonElastic(){}
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element*anEle)
   {
     G4bool result = false;
     if(( 14 == anEle->GetZ()) &&
        ( aPart->GetKineticEnergy()<150*MeV) &&
	( G4Proton::Proton() == aPart->GetDefinition() ) ) result = true;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle* aPart, const G4Element*, G4double aTemp)
   {
      G4double result = 0;
      const G4double kineticEnergy = aPart->GetKineticEnergy()/eV;
      result += the28Abundance*the28Data.getCrossSection(kineticEnergy);
      result += the29Abundance*the29Data.getCrossSection(kineticEnergy);
      result += the30Abundance*the30Data.getCrossSection(kineticEnergy);
      return result;
   }

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "CCalSiliconNonElastic: uses ADL data"<<endl;}
   
   private:
   // abundances are in atom%
   // energiesin eV; cross-sections in barns
   
   const double the28Abundance;
   CCalDataSet  the28Data;
   const double the29Abundance;
   CCalDataSet  the29Data;
   const double the30Abundance;
   CCalDataSet  the30Data;

};


#endif
