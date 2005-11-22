//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#ifndef G4NeutronHPorLFissionModel_h
#define G4NeutronHPorLFissionModel_h 1

#include "G4HadronicInteraction.hh"
#include "G4NeutronHPorLFission.hh"
#include "G4NeutronHPNames.hh"
#include "G4LFission.hh"

class G4NeutronHPorLFissionModel : public G4HadronicInteraction
{
   public:
      G4NeutronHPorLFissionModel();
      ~G4NeutronHPorLFissionModel();

      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);
      G4VCrossSectionDataSet* GiveHPXSectionDataSet() { return theHPFission->GiveXSectionDataSet(); } 

   private: 
      G4NeutronHPorLFission* theHPFission;
      G4LFission* theLFission;

      G4NeutronHPNames* theHPNames; 
};
#endif
