#ifndef CMSHydrogenNonElastic_h
#define CMSHydrogenNonElastic_h

// by JPW, working, but to be cleaned up. @@@@

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4Proton.hh"

#include "CMSDataSet.hh"

class CMSHydrogenNonElastic : public G4VCrossSectionDataSet
{
   public:
   
   CMSHydrogenNonElastic();
   virtual ~CMSHydrogenNonElastic() {}
   
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
      result = theH2CaptureData.getCrossSection(kineticEnergy);
      result += theH2DissociationData.getCrossSection(kineticEnergy);
      result *=theH2Abundance;
      return result;
   }

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "CMSHydrogenNonElastic: uses ADL data"<<endl;}
   
   private:
   // abundances are in atom%
   // energies in eV; cross-sections in barns
   
   const double theH2Abundance;
   CMSDataSet theH2CaptureData;
   CMSDataSet theH2DissociationData;

};

#endif
