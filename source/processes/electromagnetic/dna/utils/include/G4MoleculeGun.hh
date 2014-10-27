// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef MOLECULEGUN_HH_
#define MOLECULEGUN_HH_

#include "G4ITGun.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4Track;
class G4MoleculeGunMessenger;

class G4MoleculeGun : public G4ITGun
{
public:
	G4MoleculeGun();
	virtual ~G4MoleculeGun();

	virtual void DefineTracks();
	void AddMolecule(const G4String& name, const G4ThreeVector& position, double time = 0);
	void AddNMolecules(size_t n, const G4String& name, const G4ThreeVector& position, double time = 0);

protected:
	G4Track* BuildTrack(const G4String& name, const G4ThreeVector& position, double time = 0);
	std::vector<G4Track*> fTracks;
	G4MoleculeGunMessenger* fpMessenger;

};

#endif /* MOLECULEGUN_HH_ */
