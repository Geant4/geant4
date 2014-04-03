// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source, article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, in addition to the general paper on Geant4-DNA:
//
// The Geant4-DNA project, S. Incerti et al., Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following papers reference papers on chemistry:
//
// Diﬀusion-controlled reactions modelling in Geant4-DNA, M. Karamitros et al., 2014 (submitted)
// Modeling Radiation Chemistry in the Geant4 Toolkit, M. Karamitros et al., Prog. Nucl. Sci. Tec. 2 (2011) 503-508

#ifndef G4VUSERCHEMISTRYLIST_HH_
#define G4VUSERCHEMISTRYLIST_HH_

class G4Molecule;
class G4DNAMolecularReactionTable;
class G4VITModel;
class G4MoleculeDefinition;

class G4VUserChemistryList {
public:
	G4VUserChemistryList();
	virtual ~G4VUserChemistryList();

	////////////////////////////////
	// to be called from PhysicsList

	virtual void ConstructMolecule(){;} // PhysicsList::ConstructParticle
	virtual void ConstructProcess(){;}  // PhysicsList::ConstructProcess

	/////////////

	virtual void ConstructDissociationChannels() {;}
	virtual void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable) = 0;
	virtual void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable) = 0;

	void BuildPhysicsTable();

protected:
	void RegisterTimeStepModel(G4VITModel* timeStepModel, double startingTime = 0);
	void BuildPhysicsTable(G4MoleculeDefinition*);

	int verboseLevel;
};

#endif /* G4VUSERCHEMISTRYLIST_HH_ */
