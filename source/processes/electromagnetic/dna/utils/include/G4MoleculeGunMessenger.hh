// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source, article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following papers reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef MOLECULEGUNMESSENGER_HH_
#define MOLECULEGUNMESSENGER_HH_

#include "G4UImessenger.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4MoleculeGun;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIdirectory;

class G4MoleculeGunMessenger : public G4UImessenger {
public:
	G4MoleculeGunMessenger();
	virtual ~G4MoleculeGunMessenger();

	virtual void SetNewValue(G4UIcommand * command,G4String newValue);
	virtual G4String GetCurrentValue(G4UIcommand * command);
	void DefineTracks(G4MoleculeGun*);

protected:
	G4UIdirectory* fpGunDir;
	G4UIcmdWithAString* fpGunNewGunType;

	struct MultipleGun : public G4UImessenger
	{
		MultipleGun(const G4String& name, G4MoleculeGunMessenger*);
		virtual ~MultipleGun();
		virtual void SetNewValue(G4UIcommand * command,G4String newValue);
		virtual G4String GetCurrentValue(G4UIcommand * command);
		void DefineTracks(G4MoleculeGun*);

		G4UIdirectory* fpGunType;
		G4UIcmdWithAString* fpGunMoleculeModel;
		G4UIcmdWith3VectorAndUnit* fpGunPosition;
		G4UIcmdWithADoubleAndUnit* fpGunTime;
		G4UIcmdWithAnInteger* fpGunN;

		G4String fMoleculeName;
		G4ThreeVector fPosition;
		G4double fTime;
		G4int fNumber;
	};

	MultipleGun* CreateNewType(const G4String& name);

	std::vector<MultipleGun*> fMultipleGun;
};

#endif /* MOLECULEGUNMESSENGER_HH_ */
