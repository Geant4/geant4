#include "G4MWeightWindowConfigurator.hh"
#include "G4WeightWindowAlgorithm.hh"

G4MWeightWindowConfigurator::
G4MWeightWindowConfigurator(const G4String &particlename,
			    G4VWeightWindowStore &wwstore,
			    const G4VWeightWindowAlgorithm *wwAlg,
			    G4PlaceOfAction placeOfAction) :
  fPlacer(particlename),
  fWeightWindowStore(wwstore),
  fDeleteWWalg( ( ! wwAlg) ),
  fWWalgorithm(( (fDeleteWWalg) ? 
		 new G4WeightWindowAlgorithm(5,3,5) : wwAlg)),
  fMassWeightWindowProcess(0),
  fPlaceOfAction(placeOfAction)
{}

G4MWeightWindowConfigurator::
~G4MWeightWindowConfigurator()
{  
  if (fMassWeightWindowProcess) {
    fPlacer.RemoveProcess(fMassWeightWindowProcess);
    delete fMassWeightWindowProcess;
  }
  if (fDeleteWWalg) {
    delete fWWalgorithm;
  }
}

void G4MWeightWindowConfigurator::
Configure(G4VSamplerConfigurator *preConf)
{
  const G4VTrackTerminator *terminator = 0;
  if (preConf) {
    terminator = preConf->GetTrackTerminator();
  };

  fMassWeightWindowProcess = 
    new G4MassWeightWindowProcess(*fWWalgorithm, 
				  fWeightWindowStore, 
				  terminator,
				  fPlaceOfAction);
  fPlacer.AddProcessAsSecondDoIt(fMassWeightWindowProcess);
}

const G4VTrackTerminator *G4MWeightWindowConfigurator::
GetTrackTerminator() const 
{
  return fMassWeightWindowProcess;
}

