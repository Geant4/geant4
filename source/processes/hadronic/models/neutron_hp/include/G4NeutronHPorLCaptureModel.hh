//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#ifndef G4NeutronHPorLCaptureModel_h
#define G4NeutronHPorLCaptureModel_h 1

#include "G4HadronicInteraction.hh"
//#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPorLCapture.hh"
#include "G4NeutronHPNames.hh"
#include "G4LCapture.hh"

class G4NeutronHPorLCaptureModel : public G4HadronicInteraction
{
   public:
      G4NeutronHPorLCaptureModel();
      ~G4NeutronHPorLCaptureModel();

      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);
      G4VCrossSectionDataSet* GiveHPXSectionDataSet() { return theHPCapture->GiveXSectionDataSet(); } 

   private: 
      G4NeutronHPorLCapture* theHPCapture;
      G4LCapture* theLCapture;

      G4NeutronHPNames* theHPNames; 
};
#endif
