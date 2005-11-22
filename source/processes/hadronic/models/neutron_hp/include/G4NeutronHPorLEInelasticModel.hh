//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#ifndef G4NeutronHPorLEInelasticModel_h
#define G4NeutronHPorLEInelasticModel_h 1

#include "G4HadronicInteraction.hh"
//#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPorLEInelastic.hh"
#include "G4NeutronHPNames.hh"
#include "G4LENeutronInelastic.hh"

class G4NeutronHPorLEInelasticModel : public G4HadronicInteraction
{
   public:
      G4NeutronHPorLEInelasticModel();
      ~G4NeutronHPorLEInelasticModel();

      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);
      G4VCrossSectionDataSet* GiveHPXSectionDataSet() { return theHPInelastic->GiveXSectionDataSet(); } 

   private: 
      //G4NeutronHPElastic* theHPElastic;
      G4NeutronHPorLEInelastic* theHPInelastic;
      G4LENeutronInelastic* theLEInelastic;

      G4NeutronHPNames* theHPNames; 
};
#endif
