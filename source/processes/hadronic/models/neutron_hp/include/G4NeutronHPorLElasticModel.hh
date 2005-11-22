//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#ifndef G4NeutronHPorLElasticModel_h
#define G4NeutronHPorLElasticModel_h 1

#include "G4HadronicInteraction.hh"
//#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPorLElastic.hh"
#include "G4NeutronHPNames.hh"
#include "G4LElastic.hh"

class G4NeutronHPorLElasticModel : public G4HadronicInteraction
{
   public:
      G4NeutronHPorLElasticModel();
      ~G4NeutronHPorLElasticModel();

      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);
      G4VCrossSectionDataSet* GiveHPXSectionDataSet() { return theHPElastic->GiveXSectionDataSet(); } 

   private: 
      //G4NeutronHPElastic* theHPElastic;
      G4NeutronHPorLElastic* theHPElastic;
      G4LElastic* theLElastic;

      G4NeutronHPNames* theHPNames; 
};
#endif
