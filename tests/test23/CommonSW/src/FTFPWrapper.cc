
#include "FTFPWrapper.hh"

#include "G4TheoFSGenerator.hh"

#include "G4FTFModel.hh"

#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"

#include "G4SystemOfUnits.hh"

FTFPWrapper::FTFPWrapper( const G4String& name, G4ProcessType type )
   : ProcessWrapper(name,type)
{

   // override
   fUseLundStrFragm = true;

}

void FTFPWrapper::Compose()
{

   G4TheoFSGenerator* gen = new G4TheoFSGenerator();

   // fInteractionModel = new G4TheoFSGenerator();

   fStringModel      = new G4FTFModel();
   
   // Note: one can explicitly instantiate it with G4PreCompoundModel,
   //       but if not, the ctor will either fish it out of a registry,
   //       or will create a new one on the spot
   //
   fCascade          = new G4GeneratorPrecompoundInterface();
   fCascade->SetCaptureThreshold(10.*MeV);
   
   // there're these 2 options for modeling string fragmentation,
   // the Lund one is "traditionally" used with FTF, and 
   // the QGSM one (older code) is typically used with QGS
   //
   if ( fUseLundStrFragm )
   {
      fStringDecay      = new G4ExcitedStringDecay(new G4LundStringFragmentation());
   }
   else
   {
      fStringDecay      = new G4ExcitedStringDecay( new G4QGSMFragmentation() );
   }
   fStringModel->SetFragmentationModel( fStringDecay );
   
   gen->SetTransport( fCascade );
   gen->SetHighEnergyGenerator( fStringModel );
   
   fInteractionModel = gen;
   
   fInteractionModel->SetMinEnergy(GeV); 
   
   return;  

}
